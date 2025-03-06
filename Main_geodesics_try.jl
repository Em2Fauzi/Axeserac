ENV["MPLBACKEND"]="Qt5Agg"
using PyPlot
pygui(:qt5)
using Roots
using Distributed
# addprocs(16)



function main_geodesics(x0,u0,Eps,t_end;dt=0.1,int_dir=1)
    # x0 = [0.0, 500.0, pi/2, 0.0]
    # u0 = [0,-1,0,-0.5]

    nt = Int64(t_end/dt)

    e = [Eps, int_dir*dt, 1e-5]
    x = x0
    u = u_cov(x,u0,e[1])

    xup = zeros(4,nt)
    xup[:,1] = x
    for i in range(2,nt)
        (x,u,) = integr(x,u,e)

        xup[:,i] = x
    end

    return xup
end


function main_raytracing(;resy = 500, resx = 500, FoV = 12*Mass(),th0=pi/2,prograde=true,dopshift=true,allshift=true,thickacc=false,r_guess=3*Mass(),err_low=1e-8,err_high=1e-6)
    Resolution          = [resy, resx];
    
    ScreenDistance      = 0.001;
    r0                  = 500;

    Proret              = 1*prograde - 1*!prograde;
    doppler_redshift    = dopshift;
    all_redshift        = allshift;
    thick_accretion     = thickacc;

    yg                  = -2;
    mu                  = disk_isco(r_guess,Proret);

    sig                 = Mass()/4;

    lower_R             = 0.9*mu;
    upper_R             = 12*Mass();

    dxdu                = 1e-6;
    t_end               = 10000000;

    # FoV                 = 12*Mass();
    # th0                 = pi/2-0.8;
    # err_low             = 1e-8;
    # err_high            = 1e-6;

    dt = -1; 
    e = [0, dt, dxdu];
    eh = [0, dt/2, dxdu];

    dr0 = 1;
    x0 = [0, r0, th0, 0];

    param = [yg, mu, sig];

    ScreenHeight_y  = ScreenDistance * FoV / r0;
    ScreenHeight_z  = ScreenDistance * FoV / r0;
    ymage           = LinRange(-ScreenHeight_y, ScreenHeight_y, Resolution[1]);
    zmage           = LinRange(-ScreenHeight_z, ScreenHeight_z, Resolution[2]);

    y_ob = sqrt(x0[2]^2 + Spin()^2)*sin(x0[4])*sin(x0[3]);
    x_ob = sqrt(x0[2]^2 + Spin()^2)*cos(x0[4])*sin(x0[3]);
    z_ob = x0[2]*cos(x0[3]);

    dz_init         =@. dr0 * zmage / ScreenDistance;
    dy_init         =@. dr0 * ymage / ScreenDistance;
    dph0            =@. dy_init / (r0*sin(x0[3]));
    dth0            =@. dz_init / r0;

    Image           = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    NegativeR       = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    Phirec          = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    Trec            = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    nt              = Int64(t_end/dt);
    @time begin
    @sync @distributed for i in range(1,Resolution[1])
        println(i)
        for j in range(1,Resolution[2])
            dt = -1;
            x = x0;

            xcor = sqrt(x[2]^2 + Spin()^2)*cos(x[4])*sin(x[3]);
            ycor = sqrt(x[2]^2 + Spin()^2)*sin(x[4])*sin(x[3]);
            zcor = x[2]*cos(x[3]);

            u_con = [dt, dr0, -dth0[j], -dph0[i]];
            u = u_cov(x,u_con,0);

            H_init = H(x,u,e);
            H_e = H_init;
            Inten = 0;

            t_iter = 0;
            iter = 1;

            breaking = false;

            phi_rec = 0;
            t_rec = 0;
            repet = 0;
            neg_r = 0;
            while t_iter<t_end
                
                (x1,u1,xh2,uh2) = integr(x,u,[0, dt, dxdu]);

                # (x1,u1) = integr(x,u,[0, dt, dxdu]);
                # (xh1,uh1) = integr(x,u,[0, dt/2, dxdu]);
                # (xh2,uh2) = integr(xh1,uh1,[0, dt/2, dxdu]);

                if xh2[2]-dxdu<=1e-3
                    neg_r = 1
                    break
                end
                err = abs.((xh2 .- x1) ./ maximum(xh2))
                err_dom = maximum(err);

                if err_dom > err_high
                    dt = dt/2;
                    repet = repet+1
                else
                    xcor_up = sqrt(xh2[2]^2 + Spin()^2)*cos(xh2[4])*sin(x[3]);
                    ycor_up = sqrt(xh2[2]^2 + Spin()^2)*sin(xh2[4])*sin(x[3]);
                    zcor_up = xh2[2]*cos(xh2[3]);

                    if zcor_up*zcor < 0
                        r_det = x[2] - (zcor * ( (xh2[2] - x[2]) / (zcor_up - zcor) ) );
                        th_det = pi/2;

                        if r_det>lower_R && r_det<upper_R
                            iplus = disk_intensity(r_det,param);

                            gtt = g_uv([0,r_det,th_det,0])[1,1];
                            gtp = g_uv([0,r_det,th_det,0])[1,4];
                            gpp = g_uv([0,r_det,th_det,0])[4,4];
                            omeg = disk_om(r_det,Proret);

                            uphi = doppler_redshift*omeg/(sqrt(-gtt -( (2*gtp + gpp*omeg)*omeg ) ));
                            ut = (gtp*uphi + sqrt(-gtt + (gtp^2 - gtt*gpp)*uphi^2))/(-gtt);

                            g = (all_redshift*1/(ut - uh2[4]*uphi/H_init)) + (!(all_redshift)*1);
                            if imag(g)!=0
                                g=0
                            end
                            Inten = real(Inten + iplus*g^4);
                            breaking = thick_accretion;

                            #-------[For further analysis]-------#
                            t_rec = t_iter;
                            H_e = H(xh2,uh2,e);
                            phi_rec = sqrt(xh2[4]^2 + xh2[3]^2);
                            #-------[For further analysis]-------#
                        end
                    end
                    xcor = xcor_up;
                    ycor = ycor_up;
                    zcor = zcor_up;

                    x = xh2; u = uh2;

                    iter = iter+1;
                    t_iter = t_iter + abs(dt);
                    
                    if 1/g_uv(x)[2,2]<=0 || x[2]>r0 || abs(x[4]) > 40*pi || x[2]<=0
                        break
                    end

                    if err_dom > err_low
                        dt = dt/2;
                    else
                        dt = dt*2;
                    end
                end

                #-------[For further analysis]-------#
                # H_err(j,i) = abs((H_e - H_init)/H_init);
                # Phirec(j,i) = phi_rec;
                # trec(j,i) = t_rec;
                #-------[For further analysis]-------#

                if breaking
                    break
                end
                
            end
            Image[j,i] = Inten;
            Phirec[j,i] = phi_rec;
            Trec[j,i] = t_rec;
            NegativeR[j,i] = neg_r;

        end
    end
    end
    return Image, Phirec, Trec, NegativeR
end

function timedelayed_raytracing(rh_c,r_isco,t_end,TimeList;resy = 500, resx = 500, FoV = 12*Mass(),th0=pi/2,prograde=true,dopshift=true,allshift=true,thickacc=false,err_low=1e-8,err_high=1e-6)
    Resolution          = [resy, resx];
    
    ScreenDistance      = 0.001*Mass();
    r0                  = 500*Mass();

    Proret              = 1*prograde - 1*!prograde;
    doppler_redshift    = dopshift;
    all_redshift        = allshift;
    thick_accretion     = thickacc;

    yg                  = -2;
    mu                  = r_isco;

    sig                 = Mass()/4;

    lower_R             = 0.9*mu;
    upper_R             = 12*Mass();

    dxdu                = 1e-6;
    

    # FoV                 = 12*Mass();
    # th0                 = pi/2-0.8;
    # err_low             = 1e-8;
    # err_high            = 1e-6;

    dt = -1;
    
    e = [0, dt, dxdu];
    eh = [0, dt/2, dxdu];

    dr0 = 1;
    x0 = [0, r0, th0, 0];

    param = [yg, mu, sig];

    ScreenHeight_y  = ScreenDistance * FoV / r0;
    ScreenHeight_z  = ScreenDistance * FoV / r0;
    ymage           = LinRange(-ScreenHeight_y, ScreenHeight_y, Resolution[1]);
    zmage           = LinRange(-ScreenHeight_z, ScreenHeight_z, Resolution[2]);

    y_ob = sqrt(x0[2]^2 + Spin()^2)*sin(x0[4])*sin(x0[3]);
    x_ob = sqrt(x0[2]^2 + Spin()^2)*cos(x0[4])*sin(x0[3]);
    z_ob = x0[2]*cos(x0[3]);

    dz_init         =@. dr0 * zmage / ScreenDistance;
    dy_init         =@. dr0 * ymage / ScreenDistance;
    dph0            =@. dy_init / (r0*sin(x0[3]));
    dth0            =@. dz_init / r0;

    Image           = SharedArray{Float64,3}(length(TimeList),Resolution[2],Resolution[1]);
    NegativeR       = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    Phirec          = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    Trec            = SharedArray{Float64,2}(Resolution[2],Resolution[1]);
    nt              = Int64(t_end/dt);
    @time begin
    @sync @distributed for i in range(1,Resolution[1])
        println(i)
        for j in range(1,Resolution[2])
            dt = -1;
            x = x0;

            xcor = sqrt(x[2]^2 + Spin()^2)*cos(x[4])*sin(x[3]);
            ycor = sqrt(x[2]^2 + Spin()^2)*sin(x[4])*sin(x[3]);
            zcor = x[2]*cos(x[3]);

            u_con = [dt, dr0, -dth0[j], -dph0[i]];
            u = u_cov(x,u_con,0);

            H_init = H(x,u,e);
            H_e = H_init;
            Inten = 0;
            Inten0 = 0;
            detect_i0 = 1;
            hor_cross_in = 0;

            t_iter = 0;
            t_cross = 0;
            time_rec = 2
            iter = 1;

            breaking = false;

            phi_rec = 0;
            repet = 0;
            neg_r = 0;
            while t_iter<t_end
                
                (x1,u1,xh2,uh2) = integr(x,u,[0, dt, dxdu]);

                # (x1,u1) = integr(x,u,[0, dt, dxdu]);
                # (xh1,uh1) = integr(x,u,[0, dt/2, dxdu]);
                # (xh2,uh2) = integr(xh1,uh1,[0, dt/2, dxdu]);

                if xh2[2]<=1e-3
                    neg_r = 1
                    break
                end
                err = abs.((xh2 .- x1) ./ maximum(xh2))
                err_dom = maximum(err);

                if err_dom > err_high
                    dt = dt/2;
                    repet = repet+1
                else
                    xcor_up = sqrt(xh2[2]^2 + Spin()^2)*cos(xh2[4])*sin(x[3]);
                    ycor_up = sqrt(xh2[2]^2 + Spin()^2)*sin(xh2[4])*sin(x[3]);
                    zcor_up = xh2[2]*cos(xh2[3]);

                    if xh2[2] < rh_c
                        detect_i0 = 0
                        hor_cross_in = 1
                        t_cross = 0
                    end

                    if xh2[2] > rh_c && hor_cross_in ==1
                        t_cross = t_iter
                        hor_cross_in = 0
                    end

                    if zcor_up*zcor < 0
                        r_det = x[2] - (zcor * ( (xh2[2] - x[2]) / (zcor_up - zcor) ) );
                        th_det = pi/2;

                        if r_det>lower_R && r_det<upper_R
                            iplus = disk_intensity(r_det,param);

                            gtt = g_uv([0,r_det,th_det,0])[1,1];
                            gtp = g_uv([0,r_det,th_det,0])[1,4];
                            gpp = g_uv([0,r_det,th_det,0])[4,4];
                            omeg = disk_om(r_det,Proret);

                            uphi = doppler_redshift*omeg/(sqrt(-gtt -( (2*gtp + gpp*omeg)*omeg ) ));
                            ut = (gtp*uphi + sqrt(-gtt + (gtp^2 - gtt*gpp)*uphi^2))/(-gtt);

                            g = (all_redshift*1/(ut - uh2[4]*uphi/H_init)) + (!(all_redshift)*1);
                            if imag(g)!=0
                                g=0
                            end

                            if detect_i0==1
                                Inten0 = Inten0 + iplus*g^4
                            else
                                Inten = Inten + iplus*g^4
                            end

                            Image[time_rec:length(TimeList),j,i] .= Inten;
                            breaking = thick_accretion;
                        end
                    end

                    time_forward = t_iter - t_cross
                    if time_forward > TimeList[time_rec]
                        time_rec = time_rec+1
                    end

                    xcor = xcor_up;
                    ycor = ycor_up;
                    zcor = zcor_up;

                    x = xh2; u = uh2;

                    iter = iter+1;
                    t_iter = t_iter + abs(dt);
                    
                    if 1/g_uv(x)[2,2]<=0 || x[2]>r0 || abs(x[4]) > 40*pi || x[2]<=0
                        break
                    end

                    if err_dom > err_low
                        dt = dt/2;
                    else
                        dt = dt*2;
                    end
                end

                #-------[For further analysis]-------#
                # H_err(j,i) = abs((H_e - H_init)/H_init);
                # Phirec(j,i) = phi_rec;
                # trec(j,i) = t_rec;
                #-------[For further analysis]-------#

                if breaking
                    break
                end
                
            end

            Image[1,j,i] = Inten0;
            Image[2:length(TimeList),j,i] = Image[2:length(TimeList),j,i] .+ Inten0;
            NegativeR[j,i] = neg_r;

        end
    end
    end
    return Image, NegativeR
end
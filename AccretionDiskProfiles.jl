#=
Accretion disk profile
=#

#Omega
function disk_om(r,proret)
    return (-dgtp(r) + proret*sqrt(Complex(dgtp(r)^2 - dgpp(r)*dgtt(r))))/(dgpp(r));
end

#Intensity Profile
function disk_intensity(r,param)

    yg = param[1]; mu = param[2]; sig = param[3];
    
    return exp(-((yg + asinh((r-mu)/sig)).^2)/2)./(sqrt(((r-mu).^2)+sig^2));
end

#Finding R_ISCO
function disk_isco(r_guess,proret)
    L(r) = (g_uv([0,r,pi/2,0])[4,4]*disk_om(r,proret) + g_uv([0,r,pi/2,0])[1,4]);
    E(r) = -(g_uv([0,r,pi/2,0])[1,4]*disk_om(r,proret) + g_uv([0,r,pi/2,0])[1,1]);

    find_isco(r) = real(ddgtt(r)*L(r)^2 + 2*ddgtp(r)*L(r)*E(r) + ddgpp(r) * E(r)^2 - 2*(-g_uv([0,r,pi/2,0])[4,4]*disk_om(r,proret)^2 - 2*g_uv([0,r,pi/2,0])[1,4]*disk_om(r,proret) - g_uv([0,r,pi/2,0])[1,1]));

    rr = LinRange(0,10,1000)
    fisc = zeros(length(rr))

    for i in range(1,length(rr))
        fisc[i] = find_isco(rr[i])
    end
    condi = true

    rg = r_guess
    r_isco = find_zero(find_isco,rg);

    println(r_isco)
    return r_isco
end

function check_isco(r_guess,proret; limx=[0,10],limy=[-1,1])
    L(r) = (g_uv([0,r,pi/2,0])[4,4]*disk_om(r,proret) + g_uv([0,r,pi/2,0])[1,4]);
    E(r) = -(g_uv([0,r,pi/2,0])[1,4]*disk_om(r,proret) + g_uv([0,r,pi/2,0])[1,1]);

    find_isco(r) = real(ddgtt(r)*L(r)^2 + 2*ddgtp(r)*L(r)*E(r) + ddgpp(r) * E(r)^2 - 2*(-g_uv([0,r,pi/2,0])[4,4]*disk_om(r,proret)^2 - 2*g_uv([0,r,pi/2,0])[1,4]*disk_om(r,proret) - g_uv([0,r,pi/2,0])[1,1]));

    rr = LinRange(0,10,1000)
    fisc = zeros(length(rr))

    for i in range(1,length(rr))
        fisc[i] = find_isco(rr[i])
    end
    condi = true

    rg = r_guess
    r_isco = find_zero(find_isco,rg);

    fig, ax = subplots()
    ax.plot(rr,fisc)
    ax.scatter(r_isco,find_isco(r_isco))
    ax.set_xlim(limx[1],limx[2])
    ax.set_ylim(limy[1],limy[2])
    ax.grid()

    println("r_ISCO = ",r_isco)
end


#Derivative of the metric wrt r
function dgpp(r)
    dr = 1e-4;
    return (g_uv([0,r+dr,pi/2,0])[4,4]-g_uv([0,r-dr,pi/2,0])[4,4])/(2*dr);
end

function dgtt(r)
    dr = 1e-4;
    return (g_uv([0,r+dr,pi/2,0])[1,1]-g_uv([0,r-dr,pi/2,0])[1,1])/(2*dr);
end

function dgtp(r)
    dr = 1e-4;
    return (g_uv([0,r+dr,pi/2,0])[1,4]-g_uv([0,r-dr,pi/2,0])[1,4])/(2*dr);
end

function ddgpp(r)
    dr = 1e-4;
    return (dgpp(r+dr)-dgpp(r-dr))/(2*dr);
end

function ddgtt(r)
    dr = 1e-4;
    return (dgtt(r+dr)-dgtt(r-dr))/(2*dr);
end

function ddgtp(r)
    dr = 1e-4;
    return (dgtp(r+dr)-dgtp(r-dr))/(2*dr);
end


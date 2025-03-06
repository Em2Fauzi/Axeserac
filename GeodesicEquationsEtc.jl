#=

The geodesic equations and other transformation

=#


#----------------------------------------------#
# Hamiltonian related
#----------------------------------------------#
#Hamiltonian
function H(x,u,e)
    ut = u[1]; ur = u[2]; uth = u[3]; uph = u[4];

    ep = e[1];
    g = g_uv(x);

    alph =-g[1,1] + (g[1,4]^2) / (g[4,4] );
    yuu =((ur^2)/g[2,2]) + ((uth^2)/g[3,3]) + ((uph^2)/g[4,4]);

    betu =g[1,4]*uph/g[4,4];

    H=1
    if alph*(yuu + ep)<0
        H=0
        alph=0
    end

    Hamilton =H*(( sqrt(alph*(yuu + ep)) ) - betu);
    return Hamilton
end

#Hamiltonian derivative wrt x
function dHdx(x,u,e)
    dx = e[3];

    dxr = [x[1], x[2] + dx, x[3], x[4]];
    dxth = [x[1], x[2], x[3] + dx, x[4]];
    dxph = [x[1], x[2], x[3], x[4] + dx];

    d = [0,
        (H(dxr,u,e) - H(x,u,e))/dx,
        (H(dxth,u,e) - H(x,u,e))/dx,
        (H(dxph,u,e) - H(x,u,e))/dx];
    return d
end

#Hamiltonian derivative wrt u
function dHdu(x,u,e)
    du = e[3];
    
    dur = [u[1], u[2] + du, u[3], u[4]];
    duth = [u[1], u[2], u[3] + du, u[4]];
    duph = [u[1], u[2], u[3], u[4] + du];
    
    d = [1,
        (H(x,dur,e) - H(x,u,e))/du,
        (H(x,duth,e) - H(x,u,e))/du,
        (H(x,duph,e) - H(x,u,e))/du];
    return d
end

#----------------------------------------------#
# Covariant transformation
#----------------------------------------------#
function u_cov(x,u,ep)

    guv0 = g_uv(x);
    
    alph = sqrt(-guv0[1,1] + (guv0[1,4]^2) / (guv0[4,4] ));
    
    ur = guv0[2,2] * u[2];
    uth = guv0[3,3] * u[3];
    
    uph1 = -alph^2*guv0[4,4]^2*u[4];
    uph2 = guv0[4,4]*(guv0[1,4]^2 - alph^2*guv0[4,4]);
    uph3 =  (-guv0[1,4]^2*((ur^2/guv0[2,2]) + 
            (uth^2/guv0[3,3]) + ep)) + 
            alph^2*guv0[4,4]^2*u[4]^2;
    uph4 = alph^4*guv0[4,4]^4*u[4]^2;
    uph5 = guv0[1,4]^2 - alph^2*guv0[4,4];
    uph = (uph1 + sqrt(abs(uph2*uph3 + uph4)))/uph5;
    
    up = [0, ur, uth, uph];
    return up
end

#----------------------------------------------#
# Spherically symmetric
#----------------------------------------------#

function Vr(x)
    return -g_uv(x)[1,1]/g_uv(x)[3,3]
end

function ddx(x,b)
    return 1 ./ (g_uv(x)[3,3]*sqrt(( (1/b^2) - Vr(x))/(-g_uv(x)[1,1]*g_uv(x)[2,2])))
end
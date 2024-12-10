#= Integrator

Fourth-order Runge-Kutta

=#


#Integrators for Axially Symmetric
function RK4(x,u,e)
    dt = e[2];
            
    k1 =dt.*dHdu(x,u,e);
    l1 =-dt.*dHdx(x,u,e);
    
    k2 =dt.*dHdu(x+k1/2,u+l1/2,e);
    l2 =-dt.*dHdx(x+k1/2,u+l1/2,e);
    
    k3 =dt.*dHdu(x+k2/2,u+l2/2,e);
    l3 =-dt.*dHdx(x+k2/2,u+l2/2,e);
    
    k4 =dt.*dHdu(x+k3,u+l3,e);
    l4 =-dt.*dHdx(x+k3,u+l3,e);
    
    xpone = x .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4)./6;
    upone = u .+ (l1 .+ 2*l2 .+ 2*l3 .+ l4)./6;

    return (xpone,upone)
end

function RKF45(x,u,e)
    dt = e[2];
            
    k1 =dt.*dHdu(x,u,e);
    l1 =-dt.*dHdx(x,u,e);
    
    k2 =dt.*dHdu(x .+ (k1*1/4),u .+(l1*1/4),e);
    l2 =-dt.*dHdx(x .+ (k1*1/4),u .+(l1*1/4),e);

    k3 =dt.*dHdu(x .+ (k1*3/32 .+ k2*9/32),u .+(l1*3/32 .+ l2*9/32),e);
    l3 =-dt.*dHdx(x .+ (k1*3/32 .+ k2*9/32),u .+(l1*3/32 .+ l2*9/32),e);

    k4 =dt.*dHdu(x .+ (k1*1932/2197 .- k2*7200/2197 .+ k3*7296/2197),u .+(l1*1932/2197 .- l2*7200/2197 .+ l3*7296/2197),e);
    l4 =-dt.*dHdx(x .+ (k1*1932/2197 .- k2*7200/2197 .+ k3*7296/2197),u .+(l1*1932/2197 .- l2*7200/2197 .+ l3*7296/2197),e);

    k5 =dt.*dHdu(x .+ (k1*439/216 .- k2*8 .+ k3*3680/513 .- k4*845/4104),u .+(l1*439/216 .- l2*8 .+ l3*3680/513 .- l4*845/4104),e);
    l5 =-dt.*dHdx(x .+ (k1*439/216 .- k2*8 .+ k3*3680/513 .- k4*845/4104),u .+(l1*439/216 .- l2*8 .+ l3*3680/513 .- l4*845/4104),e);

    k6 =dt.*dHdu(x .+ (-k1*8/27 .+ k2*2 .- k3*3544/2565 .+ k4*1859/4104 .- k5*11/40),u .+(-l1*8/27 .+ l2*2 .- l3*3544/2565 .+ l4*1859/4104 .- l5*11/40),e);
    l6 =-dt.*dHdx(x .+ (-k1*8/27 .+ k2*2 .- k3*3544/2565 .+ k4*1859/4104 .- k5*11/40),u .+(-l1*8/27 .+ l2*2 .- l3*3544/2565 .+ l4*1859/4104 .- l5*11/40),e);
    
    
    xpone = x .+ (k1*25/216 .+ k3*1408/2565 .+ k4*2197/4104 .- k5*1/5);
    upone = u .+ (l1*25/216 .+ l3*1408/2565 .+ l4*2197/4104 .- l5*1/5);

    xcor = x .+ (k1*16/135 .+ k3*6656/12825 .+ k4*28561/56430 .- k5*9/50 .+ k6*2/55);
    ucor = u .+ (l1*16/135 .+ l3*6656/12825 .+ l4*28561/56430 .- l5*9/50 .+ l6*2/55);

    return (xpone,upone,xcor,ucor)
end


#Integrators for Spherically Symmetric

function RK4_SS(x,dx,bt)
    k1 = dx.*ddx([0,x[2],0,0],bt[1]);
    k2 = dx.*ddx([0,x[2]+(bt[2]*dx/2),0,0],bt[1]);
    k3 = dx.*ddx([0,x[2]+(bt[2]*dx/2),0,0],bt[1]);
    k4 = dx.*ddx([0,x[2]+(bt[2]*dx),0,0],bt[1]);

    return x[4] + (k1+2*k2+2*k3+k4)/6
end

function RKF45_SS(x,dx,bt)
    k1 = dx.*ddx([0,x[2],0,0],bt[1]);
    k2 = dx.*ddx([0,x[2]+(bt[2]*dx/2),0,0],bt[1]);
    k3 = dx.*ddx([0,x[2]+(bt[2]*dx/2),0,0],bt[1]);
    k4 = dx.*ddx([0,x[2]+(bt[2]*dx),0,0],bt[1]);

    return x[4] + (k1+2*k2+2*k3+k4)/6
end
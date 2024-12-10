#= Spacetime

List of spacetime

=#


#----------------------------------------------#
#Spacetime Definition
#----------------------------------------------#
function KerrNewman(x)
    M = Mass();
    a = Spin();
    Q = Charge();

    t = x[1]; r = x[2]; th = x[3]; ph = x[4];
            
    Del = r^2 + a^2 - 2*M*r + Q^2;
    Sig = r^2 + ((a^2) * (cos(th)^2));
    
    gtt = -( 1 - (2*M*r - Q^2)/Sig );
    gtph = -(a*(2*M*r - Q^2)*(sin(th)^2))/Sig;
    gphph = (r^2 + a^2 + ((2*M*r - Q^2)*(a^2)*(sin(th)^2)/Sig) )*sin(th)^2;
    
    guv =   [   [gtt,   0,          0,      gtph    ];;
                [0,     Sig/Del,    0,      0       ];;
                [0,     0,          Sig,    0       ];;
                [gtph,  0,          0,      gphph   ]];
    return guv
end

function Ghosh(x)
    M = Mass();
    a = Spin();
    Q = Charge();
    k = NewScale();

    t = x[1]; r = x[2]; th = x[3]; ph = x[4];

    mr = M*exp(-k/r)
    Del = r^2 + a^2 - 2*mr*r;
    Sig = r^2 + ((a^2) * (cos(th)^2));
    A = (r^2 + a^2)^2 - Del*a^2*sin(th)^2

    gtt = -( 1 - (2*mr*r)/Sig );
    gtph = -(a*(2*mr*r)*(sin(th)^2))/Sig;
    gphph = A*sin(th)^2/Sig
    
    guv =   [   [gtt,   0,          0,      gtph    ];;
                [0,     Sig/Del,    0,      0       ];;
                [0,     0,          Sig,    0       ];;
                [gtph,  0,          0,      gphph   ]];
    return guv
end

function Schwarzschild(x)
    M = Mass();

    t = x[1]; r = x[2]; th = x[3]; ph = x[4];
    
    gtt = -( 1 - (2*M/r) );
    
    guv =   [   [gtt,   0,      0,      0               ];;
                [0,     -1/gtt, 0,      0               ];;
                [0,     0,      r^2,    0               ];;
                [0,     0,      0,      r^2 * sin(th)^2 ]];
    return guv
end
%% Mars atmospheric drag model
function D = drag(r,v,R,S_D,C_D)   
    r_total = norms(r); %Radial position
    v_total = norms(v); %Total velocity magnitude
    P = 0.699*exp(-0.00009*(r_total-R)); %Atmospheric pressure, P = f(r)
    Temp = -31-0.000998*(r_total-R); %Atmospheric temperature, T = f(r)
    rho = P./(0.1921*(Temp+273.1)); %Atmospheric density, rho = f(P,T)
    D = 0.5*S_D*C_D*rho.*v_total; %Aerodynamic drag
end
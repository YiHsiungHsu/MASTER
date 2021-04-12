function dy = calculate_CDMTheta(u_m,Theta,A_m,f1,kappa,kappa21)
F = (0.75+0.75*(kappa+kappa21/(1+kappa21))*sqrt(u_m/A_m)...
    +2*u_m^2*f1*(1-cos(Theta(1))))/(0.75-u_m^2*f1*Theta(2)^2);
C = 1/F;

     dy = zeros(2,1); % a column vector
     dy(1) = Theta(2);
     dy(2) = -Theta(2)/u_m*(5/2-2*u_m^2*f1*Theta(2)^2-...
         (1+3/4*(kappa+(1-kappa)*kappa21/(1+kappa21))*(u_m./A_m).^(1/2))*C)...
         -sin(Theta(1))*C;
    

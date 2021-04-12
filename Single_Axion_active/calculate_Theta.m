function dy = calculate_Theta(u_m,Theta,A_m,beta_2,kappa)
C = (0.75-beta_2*u_m.^2.*Theta(2).^2)./(3/4+kappa*3/4*(u_m./A_m).^(1/2)+2*beta_2*u_m.^2.*(1-cos(Theta(1))));

     dy = zeros(2,1); % a column vector
     dy(1) = Theta(2);
     dy(2) = -Theta(2)/u_m*(5/2-2*beta_2*u_m^2*Theta(2)^2-(1+3/4*kappa*(u_m./A_m).^(1/2))*C)...
         -sin(Theta(1))*C;
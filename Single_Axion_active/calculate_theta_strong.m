function dy = calculate_theta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,beta_2,kappa,coef_photon)
%psi = [Theta, Theta', Delta_ph , g_ph, theta, theta', nu_turnon]
C = (0.75-beta_2*u_m.^2.*theta(2).^2)./(3/4+3/4*kappa*(u_m./A_m).^(1/2)...
    +2*beta_2*u_m.^2.*(1-cos(theta(1))));
F = 1/C;
if theta(7) == 1
    Phi = (u_m*F*theta(2)*theta(6)+u_m*sin(theta(1))*theta(5)+1.5*F*theta(2)*theta(5)+...
        3/8/beta_2/u_m*theta(3)*(1+3/4*kappa*sqrt(u_m/A_m)))/...
        (u_m*F*theta(2)^2-A_k/beta_2);
    Phi_p_term = beta_2*theta(2)*theta(5)-theta(4)*C/4/u_m*(1+3/4*kappa*sqrt(u_m/A_m));
    %This is the term d(phi)/d(u_m)+phi/(2u_m)
else
    Phi = (u_m*F*theta(2)*theta(6)+u_m*sin(theta(1))*theta(5)+1.5*F*theta(2)*theta(5)+...
        3/8/beta_2/u_m*theta(3)*(1/(1+alpha_nu)+3/4*kappa*sqrt(u_m/A_m)))/...
        (u_m*F*theta(2)^2-A_k/beta_2);
    Phi_p_term = beta_2*theta(2)*theta(5)-theta(4)*C/4/u_m*(1/(1+alpha_nu)+3/4*kappa*sqrt(u_m/A_m));
end
if u_m<=1/A_k/10
    u_k = sqrt(4*A_k*u_m/3);
    Phi = -coef_photon*(-cos(u_k)/u_k^2+sin(u_k)/u_k^3)/2;
end
const_Del = theta(3)/2/u_m/(1+3/4*alpha_b*sqrt(u_m/A_m));
const_g = 2*beta_2*u_m*theta(2)^2*theta(4)+C/u_m*(1+3/4*kappa*sqrt(u_m/A_m))*theta(4);

     dy = zeros(7,1); % a column vector
     dy(1) = theta(2);
     dy(2) = -theta(2)/u_m*(5/2-2*beta_2*u_m^2*theta(2)^2-(1+3/4*kappa*(u_m./A_m).^(1/2))*C)...
         -sin(theta(1))*C;
     dy(3) = 2*A_k*C/3*theta(4)+4*Phi_p_term+const_Del+const_g;
     dy(4) = -const_Del-const_g-2*Phi/u_m;
     dy(5) = theta(6);
     dy(6) = -theta(6)/u_m*(2.5-2*beta_2*u_m^2*theta(2)^2-C*(1+3/4*kappa*sqrt(u_m/A_m)))-...
         C*theta(5)*(A_k/u_m+cos(theta(1)))...
         +4*Phi_p_term*theta(2)-2*Phi*(theta(2)/u_m+C*sin(theta(1)));
     if theta(7) == 1
         if theta(3)<=0
             dy(7) = 100000;
         else
             dy(7) = 0;
         end
     else
         dy(7) = 0;
     end
function dy = calculate_theta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
    alpha_nu,epsilon0,beta_2,coef_photon,kappa)
%psi = [Delta_ph, Del_ph_p ,phic ,phis,nu_turnon,Delta_b,Del_b_p]
F = 1 + kappa*sqrt(u_m/A_m)+8/3*beta_2*sqrt(u_m)*epsilon0;
C = 1/F;
u_k = sqrt(4/3*A_k*u_m);
delta = sqrt(1-epsilon0*u_m^(-1.5)/4);
if theta(7) == 1
    Phi = -epsilon0*beta_2/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)+kappa*sqrt(u_m/A_m)*theta(6));
else
    Phi = -epsilon0*beta_2/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)/(1+alpha_nu)+kappa*sqrt(u_m/A_m)*theta(6));
end
phi_p = -coef_photon*(sin(u_k)/u_k^2+3*cos(u_k)/u_k^3-3*sin(u_k)/u_k^4)*sqrt(1/3*A_k/u_m)/2;

     dy = zeros(7,1); % a column vector
     dy(1) = theta(2);
     dy(2) = -(alpha0/(alpha0+beta0/sqrt(u_m)))/2/u_m*(A_k/3/sqrt(F)*u_m^1.5/...
         (alpha0+beta0/sqrt(u_m))+1)*theta(2)-(3-C)/4/u_m*theta(2)...
         -A_k/3*beta0/(alpha0*sqrt(u_m)+beta0)/u_m/F*theta(1)...
         -4/3*Phi*A_k/u_m*C;
     dy(3) = A_k/2/u_m/delta*sqrt(C)*theta(4)+2*phi_p/delta...
         +3/32/u_m^2/delta/sqrt(F)*theta(4);
     dy(4) = -A_k/2/u_m/delta*sqrt(C)*theta(3)...
         -Phi*sqrt(C)/delta-3/32/u_m^2/delta/sqrt(F)*theta(3)...
         +epsilon0/4/sqrt(F)/delta/u_m^1.5*theta(3);
     if theta(5) == 1
         if theta(1)<=0
             dy(5) = 1000;
         else
             dy(5) = 0;
         end
     else
         dy(5) = 0;
     end
     dy(6) = theta(7);
     dy(7) = -(3-C)*theta(7)/4/u_m-(A_k/8/sqrt(F)*beta0/(alpha0+beta0/sqrt(u_m))^2+...
         3/8/u_m*alpha0/(alpha0+beta0/sqrt(u_m)))*theta(2)-...
         A_k/4/F/u_m*beta0/(alpha0*sqrt(u_m)+beta0)*theta(1)-...
         A_k/F/u_m*Phi;


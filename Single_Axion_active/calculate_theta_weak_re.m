function dy = calculate_theta_weak_re(u_m,theta,A_k,A_m,A_t,alpha_b,alpha_nu,epsilon0,beta_2,...
    m_e,BE,T_CMB,u_mn,n_b0,kappa)
%psi = [Delta_ph, g_ph ,phic ,phis,nu_turnon,Delta_b,Q]

%calculate At after CMB
T = T_CMB*sqrt(u_mn/u_m);
n_b =  n_b0*(u_mn/u_m)^1.5;
P = (m_e*T/2/pi)^1.5*exp(-BE/T)/n_b;
X = (-P+sqrt(P^2+4*P))/2;
if abs(X) > 1
    X = 1;
end
A_t = A_t*X;
%calculate At after CMB

F = 1 + kappa*sqrt(u_m/A_m)+8/3*beta_2*sqrt(u_m)*epsilon0;
delta = sqrt(1-epsilon0*u_m^(-1.5)/4);
C = 1/F;
g_b = 3/4*(theta(2)+theta(7));
if theta(7) == 1
    Phi = -epsilon0*beta_2/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)+kappa*sqrt(u_m/A_m)*theta(6));
else
    Phi = -epsilon0*beta_2/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)/(1+alpha_nu)+kappa*sqrt(u_m/A_m)*theta(6));
end
const_Delr = theta(1)/2/u_m;
const_gr = 1/4/u_m*(3+C)*theta(2);
const_Atr = A_t/u_m^(1.5)*sqrt(C)*theta(7);
const_gb = 1/4/u_m*(3+C)*g_b;
const_Atb = A_t/alpha_b/u_m^2*sqrt(A_m*C)*theta(7);
const_gQ = 1/4/u_m*(3+C)*theta(7);
const_AtQ = (4/3*sqrt(A_m/u_m)/alpha_b+1)/u_m^1.5*A_t*sqrt(C)*theta(7);

     dy = zeros(7,1); % a column vector
     dy(1) = 2*A_k*C/3*theta(2)+const_Delr+const_gr-const_Atr;
     dy(2) = -const_Delr-const_gr-2*Phi/u_m+const_Atr;
     dy(3) = A_k/2/u_m/delta*sqrt(C)*theta(4)+3/32/u_m^2/delta/sqrt(F)*theta(4);
     dy(4) = -A_k/2/u_m/delta*sqrt(C)*theta(3)...
         -Phi*sqrt(C)/delta-3/32/u_m^2/delta/sqrt(F)*theta(3)...
         +epsilon0/4/sqrt(F)/delta/u_m^1.5*theta(3);
     if theta(5) == 1
         if theta(1)<=0
             dy(5) = 100000;
         else
             dy(5) = 0;
         end
     else
         dy(5) = 0;
     end
     dy(6) = const_gb+2*A_k*C*g_b/3+const_Atb;
     dy(7) = const_Delr-const_gQ-const_AtQ;
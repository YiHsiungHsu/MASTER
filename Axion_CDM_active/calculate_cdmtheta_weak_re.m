function dy = calculate_cdmtheta_weak_re(u_m,theta,A_k,A_m,A_t,alpha_b,alpha_nu...
    ,epsilon01,f1,m_e,BE,T_CMB,u_mn,n_b0,kappa,kappa21)
%psi = [Delta_ph, g_ph ,phic1 ,phis1, Delta_cdm,g_c,nu_turnon,Delta_b,Q]

%calculate At after CMB
T = T_CMB*sqrt(u_mn/u_m);
n_b =  n_b0*(u_mn/u_m)^1.5;
P = (m_e*T/2/pi)^1.5*exp(-BE/T)/n_b;
X = (-P+sqrt(P^2+4*P))/2;
A_t = A_t*X;
%calculate At after CMB

F = 1 + (kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m/A_m)+8/3*f1*sqrt(u_m)*epsilon01;
C = 1/F;
delta1 = 1-epsilon01*u_m^(-1.5)/4;
g_b = 3/4*(theta(2)+theta(9));
if theta(7) == 1
    Phi = -epsilon01*f1/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)+kappa*sqrt(u_m/A_m)*theta(8)...
        +(1-kappa)*kappa21/(1+kappa21)*sqrt(u_m/A_m)*theta(5));
else
    Phi = -epsilon01*f1/A_k*u_m^(-1/2)*(2*theta(3)-3*sqrt(F)/2/u_m*theta(4))...
        -3/8/A_k/u_m*(theta(1)/(1+alpha_nu)+kappa*sqrt(u_m/A_m)*theta(8)...
        +(1-kappa)*kappa21/(1+kappa21)*sqrt(u_m/A_m)*theta(5));
end
const_Delr = theta(1)/2/u_m;
const_gr = 1/4/u_m*(3+C)*theta(2);
const_Atr = A_t/u_m^(1.5)*sqrt(C)*theta(9);
const_gb = 1/4/u_m*(3+C)*g_b;
const_Atb = A_t/alpha_b/u_m^2*sqrt(A_m*C)*theta(9);
const_gQ = 1/4/u_m*(3+C)*theta(9);
const_AtQ = (4/3*sqrt(A_m/u_m)/alpha_b+1)/u_m^1.5*A_t*sqrt(C)*theta(9);
const_gc = 1/4/u_m*(3+C)*theta(6);

     dy = zeros(9,1); % a column vector
     dy(1) = 2*A_k*C/3*theta(2)+const_Delr+const_gr-const_Atr;
     dy(2) = -const_Delr-const_gr-2*Phi/u_m+const_Atr;
     dy(3) = A_k/2/u_m/delta1*sqrt(C)*theta(4)+...
         3/32/u_m^2/delta1/sqrt(F)*theta(4);
     dy(4) = -A_k/2/u_m/delta1*sqrt(C)*theta(3)...
         -Phi*sqrt(C)/delta1-3/32/u_m^2/delta1/sqrt(F)*theta(3)...
         +epsilon01/4/sqrt(F)/delta1/u_m^1.5*theta(3);
     dy(5) = const_gc+2*A_k*C*theta(6)/3;
     dy(6) = -const_gc-3*Phi/2/u_m;
     if theta(7) == 1
         if theta(1)<=0
             dy(7) = 100000;
         else
             dy(7) = 0;
         end
     else
         dy(7) = 0;
     end
     dy(8) = const_gb+2*A_k*C*g_b/3+const_Atb;
     dy(9) = const_Delr-const_gQ-const_AtQ;
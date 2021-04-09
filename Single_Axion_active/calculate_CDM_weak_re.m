function dy = calculate_CDM_weak_re(a,Delta,A_k,alpha_nu,A_t, m_e,BE,T_CMB,aeq,n_b0,kappa)
%psi = [Delta_ph, g_ph,Delta_D,g_D,nu_turnon,Delta_b,Q]
alpha_t = (1+alpha_nu);
alpha_b = alpha_t*kappa;
alpha_D = alpha_t*(1-kappa);
T = T_CMB/aeq/a;
n_b =  n_b0/(aeq*a)^3;
P = (m_e*T/2/pi)^1.5*exp(-BE/T)/n_b;
X = (-P+sqrt(P^2+4*P))/2;
A_t = A_t*X;
g_b = 3/4*(Delta(7)+Delta(2));

if Delta(5) == 1
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2*alpha_t+(alpha_b*Delta(6)+alpha_D*Delta(3))/a);
    Phi_p_term = -(Delta(2)/a*alpha_t+alpha_b*g_b+alpha_D*Delta(4))/(2*alpha_t*(a+1));
    %This is the term d(phi)/d(u_m)+phi/(2u_m)
else
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2+(alpha_b*Delta(6)+alpha_D*Delta(3))/a);
    Phi_p_term = -(Delta(2)/a+alpha_b*g_b+alpha_D*Delta(4))/(2*alpha_t*(a+1));
end

const_gr = 3/2/a*(1+1/3/(a+1))*Delta(2);
const_Atr = A_t/a^2/sqrt(1+a)*Delta(7);
%const_Atr = A_k/3/A_t/(1/a^2+4/3/a^3/alpha_b)/(1/a+1);
const_gd = 3/2/a*(1+1/3/(a+1))*Delta(4);
const_gb = 3/2/a*(1+1/3/(a+1))*g_b;
const_Atb = A_t/alpha_b/a^3/sqrt(1+a)*Delta(7);
%const_Atb = A_k/4/A_t/(1/a^2+4/3/a^3/alpha_b)/(1/a+1);
const_Q =  3/2/a*(1+1/3/(a+1))*Delta(7);
const_AtQ = A_t/a^2/sqrt(a+1)*(1+4/3/alpha_b/a)*Delta(7);


     dy = zeros(7,1); % a column vector
     dy(1) = A_k*a/3/(a+1)*Delta(2)+4*Phi_p_term+Delta(1)/a+const_gr-const_Atr;
     dy(2) = -Delta(1)/a-4*Phi/a-const_gr+const_Atr;
     dy(3) = A_k*a/3/(a+1)*Delta(4)+3*Phi_p_term+const_gd;
     dy(4) = -3*Phi/a-const_gd;
     if Delta(5) == 1
         if Delta(1)<=0
             dy(5) = 100000;
         else
             dy(5) = 0;
         end
     else
         dy(5) = 0;
     end
     dy(6) = A_k*a/3/(a+1)*g_b+3*Phi_p_term+const_gb+const_Atb;
     dy(7) = -const_Q-const_AtQ+Delta(1)/a;


function dy = calculate_CDM_weak(a,Delta,A_k,alpha_nu,A_t,kappa)
%psi = [Delta_ph, Delta_ph',Delta_D,g_D,nu_turnon,Delta_b,Delta_b']
alpha_t = (1+alpha_nu);
alpha_b = alpha_t*kappa;
alpha_D = alpha_t*(1-kappa);
%g_r = 3*(1+a)*Delta(2)/A_k/a;
%g_b = 3*(1+a)*Delta(7)/A_k/a;
a_ab = 3*a*alpha_b/(3*a*alpha_b+4);
b_ab = 4/(3*alpha_b*a+4);
if Delta(5) == 1
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2*alpha_t+(alpha_b*Delta(6)+alpha_D*Delta(3))/a);
    %Phi_p_term = -(TA*g_r/a*alpha_t+TB*alpha_b*g_b+alpha_D*Delta(4))/(2*alpha_t*(a+1));
    %This is the term d(phi)/d(u_m)+phi/(2u_m)
else
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2+(alpha_b*Delta(6)+alpha_D*Delta(3))/a);
    %Phi_p_term = -(TA*g_r/a+TB*alpha_b*g_b+alpha_D*Delta(4))/(2*alpha_t*(a+1));
end

%const_gr = 3/2/a*(1+1/3/(a+1))*g_r;
%const_Atr = A_t/a^2/sqrt(1+a)*(4/3*Delta(7)-Delta(2));
const_Atr = A_k/3/A_t/(1/a^2+4/3/a^3/alpha_b)/sqrt(a+1);
const_gd = 3/2/a*(1+1/3/(a+1))*Delta(4);
%const_gb = 3/2/a*(1+1/3/(a+1))*Delta(7);
%const_Atb = A_t/alpha_b/a^3/sqrt(1+a)*(4/3*Delta(7)-Delta(2));
const_Atb = A_k/4/A_t/(1/a^2+4/3/a^3/alpha_b)/sqrt(a+1);


     dy = zeros(7,1); % a column vector
     dy(1) = Delta(2);
     dy(2) = -a_ab*(const_Atr+1/a)*Delta(2)-1/2/(1+a)*Delta(2)...
         -A_k/3/(a+1)*b_ab*Delta(1)-4/3*A_k/(1+a)*Phi;
     dy(3) = A_k*a/3/(a+1)*Delta(4)+const_gd;
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
     dy(6) = Delta(7);
     dy(7) = -Delta(7)/2/(1+a)-const_Atb*b_ab*Delta(2)-3/4/a*a_ab*Delta(2)...
         -A_k/4*b_ab/(a+1)*Delta(1)-A_k*Phi/(1+a);


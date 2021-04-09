function dy = calculate_CDM_strong(a,Delta,A_k,alpha_nu)
%psi = [Delta_ph, g_ph,Delta_D,g_D,nu_turnon]
alpha_t = (1+alpha_nu);
alpha_b = alpha_t*0.2;
alpha_D = alpha_t*0.8;
if Delta(5) == 1
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2*alpha_t+(3/4*alpha_b*Delta(1)+alpha_D*Delta(3))/a);
    Phi_p_term = -(Delta(2)/a*alpha_t+3*alpha_b*Delta(2)/4+alpha_D*Delta(4))/(2*alpha_t*(a+1));
    %This is the term d(phi)/d(u_m)+phi/(2u_m)
else
    Phi = -1.5/A_k/alpha_t*(Delta(1)/a^2+(3/4*alpha_b*Delta(1)+alpha_D*Delta(3))/a);
    Phi_p_term = -(Delta(2)/a+3*alpha_b*Delta(2)/4+alpha_D*Delta(4))/(2*alpha_t*(a+1));
end

const_gr = 3/2/a*(1+1/3/(a+1))*Delta(2);
const_gd = 3/2/a*(1+1/3/(a+1))*Delta(4);


     dy = zeros(5,1); % a column vector
     dy(1) = A_k*a/3/(a+1)*Delta(2)+4*Phi_p_term+Delta(1)/a/(1+3*alpha_b*a/4)+const_gr;
     dy(2) = -Delta(1)/a/(1+3*alpha_b*a/4)-4*Phi/a-const_gr;
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


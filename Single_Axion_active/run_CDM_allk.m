clear all
close all

%%%input parameters
%End time in a_eq (scale factor at matter-radiation equality)
a_f = 30;

%Initial wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_initial = 0.1;

%End wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_end = 2;

%Spectral resolution in k_c22
del_k = 0.02;

%current axion background ratio
Omega_CDM = 0.24; 

%current baryon background ratio
Omega_baryon = 0.06;

%current nuetrino back ground ratio
Omega_neutrino = 3.43e-5; 

%current photon background ratio
Omega_photon = 5.04e-5; 

%current Hubble constant/100 in km/s/Mpc
Hubble_constant = 0.7;

%electron mass in eV 
m_e = 0.510998902e+6;

%neutron mass in eV 
m_n = 939.5654205254e+6;

%Thomson Scattering crossection in eV^-2
sigma_T = 1.708477740795844e-15;

%CMB current temperature
T_CMB = 2.34864e-4;

%Rydberg constant in eV
BE = 13.605693122994;

%Maximum first phase resolution in u_m
res = 0.5;
%%%input parameters

tic
%%% parameter conversion
G = 6.708830620767732e-57;
aeq = (Omega_photon+Omega_neutrino)/(Omega_CDM+Omega_baryon);
H0 = 2.133119459810834e-33*Hubble_constant;
Hr_eq = H0*sqrt((Omega_photon+Omega_neutrino)/aeq^2);
u_meq = 1e-22*aeq/2/Hr_eq;
k_c = sqrt(2*1e-22*aeq*Hr_eq);
k_c_kpc = k_c/1.97326960277e-5/3.24077929e-22;
u_mn = aeq^-2*u_meq;
rho_cr = 3*H0^2/8/pi/G;
n_b0 = rho_cr/m_n*Omega_baryon;
A_t = n_b0*sigma_T/aeq^3/...
    sqrt(8*pi*G/3*rho_cr*(Omega_baryon+Omega_CDM)/aeq^3);
k = k_initial;
n = 1;
%%% parameter conversion

%%%zeroth order parameter
alpha_nu = Omega_neutrino/Omega_photon;
kappa = Omega_baryon/(Omega_baryon+Omega_CDM);
alpha_b = (1.0+alpha_nu)*kappa;
alpha_t = (1.0+alpha_nu);
alpha_D = alpha_t*(1-kappa);
A_m = u_meq;
%%%zeroth order parameter

%%% first order parameter
coef_photon = 2/3;
u_m0 = 0.0005;
a_i = sqrt(u_m0/u_meq);
%%% first order parameter

while k<=k_end

    A_k = k^2;
    turnon_weak = 30/A_k;
    a_mid = sqrt(turnon_weak/u_meq);
    opts = odeset('MaxStep',res);
    A_kc = A_k*4*u_meq;
    a_tr = sqrt(A_t/10/sqrt(A_kc));

    %%%initial conditions
    u_k0 = sqrt(A_kc/3)*a_i;
    Del_ph0 = coef_photon*(-cos(u_k0)+sin(u_k0)/u_k0...
        -(0.25 +0.5625*alpha_b)*u_k0^2*sqrt(u_m0/A_m)/3.0); 
    g_ph0 = coef_photon*(sin(u_k0)/u_k0+2*cos(u_k0)/u_k0^2-...
        2*sin(u_k0)/u_k0^3+(1 - 0.75*alpha_b)*sqrt(u_m0/A_m)/6.0);
    %Del_ph0 = coef_photon*(2*cos(u_k0)-2*sin(u_k0)/u_k0);
    %g_ph0 = coef_photon*(-2*sin(u_k0)/u_k0-4*cos(u_k0)/u_k0^2+4*sin(u_k0)/u_k0^3);
    gamma = double(eulergamma);
    Del_D0 = 1.5*coef_photon*(sin(u_k0)/u_k0+(cos(u_k0)-1)/u_k0^2+log(u_k0)...
        -cosint(u_k0)+gamma-1/2);
    g_D0 = 1.5*coef_photon*(1/u_k0^2-sin(u_k0)/u_k0^3);
    nu_turnon = 1;
    %%%initial conditions
    
    if a_f<=a_mid
        [a_1,Delta] = ode45(@(a,Delta)...
            calculate_CDM_strong(a,Delta,A_kc,alpha_nu),...
            [a_i a_mid],[Del_ph0; g_ph0; Del_D0; g_D0; nu_turnon],opts);
        Del_b1 = 3/4.*Delta(:,1);
        g_b1 = 3/4.*Delta(:,2);
        Delta = [Delta Del_b1 g_b1];
        g_b = Delta(end,7);
        g_ph = Delta(end,2);
        
    elseif a_f<a_tr
        [a_1,Delta_1] = ode45(@(a,Delta)...
            calculate_CDM_strong(a,Delta,A_kc,alpha_nu),...
            [a_i a_mid],[Del_ph0; g_ph0; Del_D0; g_D0; nu_turnon],opts);
        Del_b1 = 3/4.*Delta_1(:,1);
        g_b1 = 3/4.*Delta_1(:,2);
        Delta_1 = [Delta_1 Del_b1 g_b1];
        
        %%%potential
        index_nuturn = find(Delta_1(:,5) == 1,1,'last');
        Delta_off = Delta_1(1:index_nuturn,:);
        a_off = a_1(1:index_nuturn,:);
        Phi_off = -1.5/A_kc/alpha_t*(Delta_off(:,1)./a_off.^2*alpha_t+...
            (3/4*alpha_b*Delta_off(:,1)+alpha_D*Delta_off(:,3))./a_off);
        Delta_on = Delta_1(index_nuturn+1:end,:);
        a_on = a_1(index_nuturn+1:end,:);
        Phi_on = -1.5/A_kc/alpha_t*(Delta_on(:,1)./a_on.^2+...
            (3/4*alpha_b*Delta_on(:,1)+alpha_D*Delta_on(:,3))./a_on);
        Phi_1 = [Phi_off;Phi_on];
        clear index_nuturn a_off Phi_off Delta_on a_on Phi_on
        %%%potential

        %%%diffusion approx IC
        diffr = gradient(Delta_1(:,1))./gradient(a_1);
        diffb = gradient(Delta_1(:,6))./gradient(a_1);
        Del_rp0 = diffr(end);
        Del_bp0 = diffb(end);
        [a_2,Delta_2] = ode45(@(a,Delta)...
            calculate_CDM_weak(a,Delta,A_kc,alpha_nu,A_t,kappa),...
            [a_mid a_f],[Delta_1(end,1); Del_rp0; Delta_1(end,3); Delta_1(end,4)...
            ;Delta_1(end,5); Delta_1(end,6); Del_bp0]);
        %%%diffusion approx IC

        if Delta_1(end,5) == 1
            index_nuturn = find(Delta_2(:,5) == 1,1,'last');
            Delta_off = Delta_2(1:index_nuturn,:);
            a_off = a_2(1:index_nuturn,:);
            Phi_off = -1.5/A_kc/alpha_t*(Delta_off(:,1)./a_off.^2*alpha_t+...
                (alpha_b*Delta_off(:,6)+alpha_D*Delta_off(:,3))./a_off);
            Delta_on = Delta_2(index_nuturn+1:end,:);
            a_on = a_2(index_nuturn+1:end,:);
            Phi_on = -1.5/A_kc/alpha_t*(Delta_on(:,1)./a_on.^2+...
                (3/4*alpha_b*Delta_on(:,6)+alpha_D*Delta_on(:,3))./a_on);
            Phi_2 = [Phi_off;Phi_on];
            clear index_nuturn a_off Phi_off Delta_on a_on Phi_on
        else
            Phi_2 = -1.5/A_kc/alpha_t*(Delta_2(:,1)./a_2.^2+...
                (3/4*alpha_b*Delta_2(:,6)+alpha_D*Delta_2(:,3))./a_2);
        end
        g_ph2 = (Delta_2(:,2)-Delta_2(:,1)./a_2)./(A_k*a_2./(a_2+1)/3+...
            3/2./a_2.*(1+1/3./(a_2+1)));
        g_b2 = Delta_2(:,7)./(A_k*a_2./(a_2+1)/3+...
            3/2./a_2.*(1+1/3./(a_2+1)));
        
        Delta = Delta_2(end,:);
        g_ph = g_ph2(end);
        g_b = g_b2(end);
        
    else
        [a_1,Delta_1] = ode45(@(a,Delta)...
            calculate_CDM_strong(a,Delta,A_kc,alpha_nu),...
            [a_i a_mid],[Del_ph0; g_ph0; Del_D0; g_D0; nu_turnon],opts);

        Del_b1 = 3/4.*Delta_1(:,1);
        g_b1 = 3/4.*Delta_1(:,2);
        Delta_1 = [Delta_1 Del_b1 g_b1];
        
        %potential
        index_nuturn = find(Delta_1(:,5) == 1,1,'last');
        Delta_off = Delta_1(1:index_nuturn,:);
        a_off = a_1(1:index_nuturn,:);
        Phi_off = -1.5/A_kc/alpha_t*(Delta_off(:,1)./a_off.^2*alpha_t+...
            (3/4*alpha_b*Delta_off(:,1)+alpha_D*Delta_off(:,3))./a_off);
        Delta_on = Delta_1(index_nuturn+1:end,:);
        a_on = a_1(index_nuturn+1:end,:);
        Phi_on = -1.5/A_kc/alpha_t*(Delta_on(:,1)./a_on.^2+...
            (3/4*alpha_b*Delta_on(:,1)+alpha_D*Delta_on(:,3))./a_on);
        Phi_1 = [Phi_off;Phi_on];
        clear index_nuturn a_off Phi_off Delta_on a_on Phi_on
        %potential

        %initial condition
        diffr = gradient(Delta_1(:,1))./gradient(a_1);
        diffb = gradient(Delta_1(:,6))./gradient(a_1);
        Del_rp0 = diffr(end);
        Del_bp0 = diffb(end);
        %initial condition
        
        [a_2,Delta_2] = ode45(@(a,Delta)...
            calculate_CDM_weak(a,Delta,A_kc,alpha_nu,A_t,kappa),...
            [a_mid a_tr],[Delta_1(end,1); Del_rp0; Delta_1(end,3);...
            Delta_1(end,4);Delta_1(end,5); Delta_1(end,6); Del_bp0]);
        
        [pkr,locr] = findpeaks(abs(Delta_2(:,1)),a_2);
        index_tr = find(abs(Delta_2(:,1))==pkr(length(locr)),1,'last');
        a_2 = a_2(1:index_tr);
        Delta_2 = Delta_2(1:index_tr,:);
        a_tr = a_2(end);

        if Delta_1(end,5) == 1
            index_nuturn = find(Delta_2(:,5) == 1,1,'last');
            Delta_off = Delta_2(1:index_nuturn,:);
            a_off = a_2(1:index_nuturn,:);
            Phi_off = -1.5/A_kc/alpha_t*(Delta_off(:,1)./a_off.^2*alpha_t+...
                (alpha_b*Delta_off(:,6)+alpha_D*Delta_off(:,3))./a_off);
            Delta_on = Delta_2(index_nuturn+1:end,:);
            a_on = a_2(index_nuturn+1:end,:);
            Phi_on = -1.5/A_kc/alpha_t*(Delta_on(:,1)./a_on.^2+...
                (3/4*alpha_b*Delta_on(:,6)+alpha_D*Delta_on(:,3))./a_on);
            Phi_2 = [Phi_off;Phi_on];
            clear index_nuturn a_off Phi_off Delta_on a_on Phi_on
        else
            Phi_2 = -1.5/A_kc/alpha_t*(Delta_2(:,1)./a_2.^2+...
                (3/4*alpha_b*Delta_2(:,6)+alpha_D*Delta_2(:,3))./a_2);
        end
        clear index_nuturn a_off Phi_off Delta_on a_on Phi_on

        %%%initial condition
        T = T_CMB/aeq/a_tr;
        n_b =  n_b0/(aeq*a_tr)^3;
        P = (m_e*T/2/pi)^1.5*exp(-BE/T)/n_b;
        X = (-P+sqrt(P^2+4*P))/2;
        A_t = A_t*X;
        
        Q_30 = sqrt(1/a_2(index_tr)^2+1/a_2(index_tr))*Delta_2(index_tr,1)/...
            (A_t/a_2(index_tr)^2+4*A_t/3/a_2(index_tr)/alpha_b+...
            sqrt(1/a_2(index_tr)^2+1/a_2(index_tr)));
        g_r30 = (Delta_2(index_tr,2)-Delta_2(index_tr,1)/a_2(index_tr)...
            +A_t/a_2(index_tr)^2/sqrt(a_2(index_tr)+1)*Q_30)...
            /(A_kc*a_2(index_tr)/3/(1+a_2(index_tr))+...
            1.5/a_2(index_tr)*(1+1/3/(a_2(index_tr)+1)));
        %%%initial condition
        
        [a_3,Delta_3] = ode45(@(a,Delta)...
            calculate_CDM_weak_re(a,Delta,A_kc,alpha_nu,A_t, m_e,BE,...
            T_CMB,aeq,n_b0,kappa),[a_tr a_f],[Delta_2(index_tr,1); ...
            g_r30; Delta_2(index_tr,3); Delta_2(index_tr,4)...
            ;Delta_2(index_tr,5); Delta_2(index_tr,6); Q_30]);

            Phi_3 = -1.5/A_kc/alpha_t*(Delta_3(:,1)./a_3.^2+...
                (alpha_b*Delta_3(:,6)+alpha_D*Delta_3(:,3))./a_3);
            g_b3 = 3/4*(Delta_3(:,7)+Delta_3(:,2));

        Delta = Delta_3(end,:);
        g_ph = Delta_3(end,2);
        g_b = g_b3(end);
    end

    kn(n) = k;
    Delta_n(n) = Delta(3);
    g_CDMn(n) = Delta(4);
    Delta_bn(n) = Delta(6);
    g_bn(n) = g_b;
    Delta_phn(n) = Delta(1);
    g_phn(n) = g_ph;

    if k == k_end
        break
    end
    
    k = k + del_k;
    if k>k_end
        k=k_end;
    end
    n = n+1;
end
toc

%save file
afstr = num2str(a_f);
ki = num2str(k_initial);
kf = num2str(k_end);
filename = ['CDM_' afstr 'aeq' '_from_k_' ki 'to' kf '.txt'];
fid = fopen(filename,'wt');
fprintf(fid,'Hubble constant now = %17.24f km/s/Mpc \n',...
    100*Hubble_constant);
fprintf(fid,'Omega_CDM = %17.24f \n',Omega_CDM);
fprintf(fid,'Omega_baryon = %17.24f \n',Omega_baryon);
fprintf(fid,'Omega_photon = %17.24f \n',Omega_photon);
fprintf(fid,'Omega_neutrino = %17.24f \n',Omega_neutrino);
fprintf(fid,'critical k_c22 = 1/(%3.5f kpc) \n',1/k_c_kpc);
fprintf(fid,'k \t CDM\t Baryon\t Photon\t g_CDM=3H theta_CDM \t g_b\t g_ph\t \n');
fprintf(fid,['%17.24f\t %17.24f\t %17.24f\t %17.24f\t %17.24f\t'...
    '%17.24f\t %17.24f\t\n'],...
    [kn;Delta_n;Delta_bn;Delta_phn;g_CDMn;g_bn;g_phn]);
fclose(fid);
%save file
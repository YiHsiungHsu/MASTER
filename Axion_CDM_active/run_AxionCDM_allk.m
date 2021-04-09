clear all
close all


%%%input parameters
%End time in a_eq (scale factor at matter-radiation equality)
a_f = 0.01;

%Initial wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_initial = 1;

%End wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_end = 1.1;

%Spectral resolution in k_c22
del_k = 0.02;

%first particle mass in 10^-22 eV
mr = 1;

%background axion angle shift from the top of the potential hill for particle 1:
angle_shift1 = 90;

%current total axion background ratio
Omega_axion = 0.24; 

%axion background ratio comp. 2 to comp. 1
kappa21 = 1;

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

%(Maximum first phase resolution)*mu in u_m
res = 0.1;

%(Maximum resolution)*mu for iteration of decay constant determination
res_i = 0.01;

%number of iteration for decay constant determination
N_i = 2;

%iteration end point in a_eq
a_stop = 0.05;

%%%input parameters

tic
%%% parameter conversion
G = 6.708830620767732e-57;
a_eq = (Omega_photon+Omega_neutrino)/(Omega_axion+Omega_baryon);
H0 = 2.133119459810834e-33*Hubble_constant;
Hr_eq = H0*sqrt((Omega_photon+Omega_neutrino)/a_eq^2);
u_meq = mr*1e-22*a_eq/2/Hr_eq;
k_c = sqrt(2*mr*1e-22*a_eq*Hr_eq);
k_c_kpc = k_c/1.97326960277e-5/3.24077929e-22;
fprintf('critical k_c = 1/(%3.5f kpc) \n',1/k_c_kpc);
u_mn = a_eq^-2*u_meq;
rho_cr = 3*H0^2/8/pi/G;
n_b0 = rho_cr/m_n*Omega_baryon;
A_t = n_b0*sigma_T*sqrt(mr*1e-22)/2/sqrt(2)/a_eq^3/...
    (8*pi*G/3*rho_cr*(Omega_baryon+Omega_axion)/a_eq^3)^0.75;
%%% parameter conversion

%%% initial setting
k = k_initial;
u_m0 = 0.0005*mr;          %starting point
u_mf = a_f^2*u_meq;    %end point
n = 1;                  % parameter for save
%%% initial setting

%%%zeroth order parameter
alpha_nu = Omega_neutrino/Omega_photon;
kappa = Omega_baryon/(Omega_baryon+Omega_axion);
alpha_b = (1.0+alpha_nu)*kappa;
A_m = u_meq;
%%%zeroth order parameter

%%% zeroth order IC
Theta01 = pi - angle_shift1/180*pi;
Theta_i1 = Theta01-sin(Theta01)*u_m0^2/5+...
    u_m0^2.5*sin(Theta01)/sqrt(A_m)*11/75;
Theta_ip1 = -2*sin(Theta01)*u_m0/5+...
    u_m0^1.5*sin(Theta01)/sqrt(A_m)*11/30;
%%% zeroth order IC

%%%calculate beta^2 by passive evolution abd Newton method
u_mstop = a_stop^2*u_meq;
[u_m_p,theta_p] = ode45(@(u_m,Theta)...
        calculate_CDMTheta_passive(u_m,Theta),...
        [u_m0 u_mstop],[Theta_i1;Theta_ip1]);
epsilon01 = (0.5*theta_p(end,2)^2+1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
%%%interation
int = 1; %number of interation
opts = odeset('MaxStep',res_i);
for i = 1:int
    [u_m_p,theta_p] = ode45(@(u_m,Theta)...
            calculate_CDMTheta(u_m,Theta,A_m,f1,kappa,kappa21),...
            [u_m0 u_mstop],[Theta_i1;Theta_ip1],opts);
    F = (3/4+3/4*(kappa+(1-kappa)*kappa21/(1+kappa21))*(u_m_p./A_m).^(1/2)...
        +2*f1*u_m_p.^2.*(1-cos(theta_p(:,1))))./...
        (0.75-f1*u_m_p.^2.*theta_p(:,2).^2);
    epsilon01 = (0.5*F(end)*theta_p(end,2)^2+...
        1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
    f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
end

clear u_m_p theta_p u_mstop epsilon0 F int i
%%%calculate beta^2 by passive evolution and Newton method

%%% first order parameter
coef_photon = 2/3;
coef_phi = coef_photon*3/2;
%%% first order parameter

%%%Thomson scattering and others
%Thomson
alpha0 = A_t/2;
beta0 = 4*sqrt(u_meq)/3/alpha_b*alpha0;
while k<=k_end 
    A_k = k^2/mr;
    u_mmid = 30/A_k;
    if u_mmid <= 50
        u_mmid = 50;
    end
    u_mtr = 2*alpha0*0.1/sqrt(A_k);

    %%%initial condition for exact strong
    u_k0 = sqrt(4/3*A_k*u_m0);
    gamma = double(eulergamma);

    Del_ph0 = coef_photon*(-cos(u_k0)+sin(u_k0)/u_k0...
        -(0.25 +0.5625*alpha_b)*u_k0^2*sqrt(u_m0/A_m)/3.0); 
    g_ph0 = coef_photon*(sin(u_k0)/u_k0+2*cos(u_k0)/u_k0^2-...
        2*sin(u_k0)/u_k0^3+(1 - 0.75*alpha_b)*sqrt(u_m0/A_m)/6.0);
    psi01 = coef_phi*sin(Theta01)*(4*u_m0^2-2*(8/15+alpha_b*21/20)...
        *u_m0^2.5/sqrt(A_m)-176*A_k*u_m0^3/210)/90;
    psi11 = coef_phi*sin(Theta01)*(8*u_m0-5*(8/15+alpha_b*21/20)*...
        u_m0^1.5/sqrt(A_m)-176*A_k*u_m0^2/70)/90;
    Del_D0 = 1.5*coef_photon*(sin(u_k0)/u_k0+(cos(u_k0)-1)/u_k0^2+log(u_k0)...
        -cosint(u_k0)+gamma-1/2);
    g_D0 = 1.5*coef_photon*(1/u_k0^2-sin(u_k0)/u_k0^3);
    nu_turnon = 1;
    %%%initial condition for exact strong

    if u_mf<=u_mmid
        %%%theta ode exact
        opts = odeset('MaxStep',res);
        [u_m,theta] = ode45(@(u_m,theta)...
            calculate_cdmtheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
            f1,kappa,kappa21),[u_m0 u_mf],...
            [Theta_i1;Theta_ip1;Del_ph0;g_ph0;...
            psi01;psi11;Del_D0;g_D0;nu_turnon],opts);
        Del_b = 3/4*theta(:,3);
        g_b = 3/4*theta(:,4);
        theta = [theta Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta(:,9) == 1,1,'last');
        F = (3/4+3/4*(kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m./A_m)+...
            2*f1*u_m.^2.*(1-cos(theta(:,1))))...
            ./(0.75-f1*u_m.^2.*theta(:,2).^2);
        theta_off = theta(1:index_nuturn,:);
        u_m_off = u_m(1:index_nuturn,:);
        F_off = F(1:index_nuturn,:);
        theta_on = theta(index_nuturn+1:end,:);
        u_m_on = u_m(index_nuturn+1:end,:);
        F_on = F(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)+...
            +3/8/f1./u_m_off.*(theta_off(:,3)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_off(:,7)...
            +kappa*theta_off(:,10)).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/f1);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)+...
            +3/8/f1./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_on(:,7)...
            +kappa*theta_on(:,10)).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/f1);
        Phi = [Phi_nuoff ; Phi_nuon];
        Phi_p = gradient(Phi)./gradient(u_m);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density and output data
        E_psi = 0.5*F(end)*theta(end,2).^2+1-cos(theta(end,1));
        del_psi = F(end)*theta(end,2)*theta(end,6)+...
            sin(theta(end,1))*theta(end,5)-Phi(end)*F(end)*theta(end,2)^2+...
            1.5*F(end)*theta(end,2)*theta(end,5)/u_m(end);
        Delta_theta1 = del_psi/E_psi;
        g_theta1 = 1.5.*F(end)*theta(end,2).*theta(end,5)/u_m(end)/E_psi;
        epsilon02 = 3/8*(1-kappa)*kappa21/(1+kappa21)/sqrt(u_meq);
        E_CDM = epsilon02*u_m(end)^-1.5;
        Delta_theta2 = theta(end,7);
        g_theta2 = theta(end,8);
        Delta_theta = (f1*E_psi.*Delta_theta1+E_CDM.*Delta_theta2)/(f1*E_psi+E_CDM);
        g_theta = (f1*E_psi*g_theta1+E_CDM*g_theta2)/(f1*E_psi+E_CDM);
        Delta_b = theta(end,10);
        g_b = theta(end,11);
        Delta_ph = theta(end,3);
        g_ph = theta(end,4);
        %%%covariant energy density and output data
    elseif u_mf<=u_mtr
        %%%theta ode exact
        opts = odeset('MaxStep',res);
        [u_m_e,theta_e] = ode45(@(u_m,theta)...
            calculate_cdmtheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
            f1,kappa,kappa21),[u_m0 u_mmid],...
            [Theta_i1;Theta_ip1;Del_ph0;g_ph0;...
            psi01;psi11;Del_D0;g_D0;nu_turnon],opts);
        Del_b = 3/4*theta_e(:,3);
        g_b = 3/4*theta_e(:,4);
        theta_e = [theta_e Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta_e(:,9) == 1,1,'last');
        F_e = (3/4+3/4*(kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m_e./A_m)+...
            2*f1*u_m_e.^2.*(1-cos(theta_e(:,1))))...
            ./(0.75-f1*u_m_e.^2.*theta_e(:,2).^2);
        theta_off = theta_e(1:index_nuturn,:);
        u_m_off = u_m_e(1:index_nuturn,:);
        F_off = F_e(1:index_nuturn,:);
        theta_on = theta_e(index_nuturn+1:end,:);
        u_m_on = u_m_e(index_nuturn+1:end,:);
        F_on = F_e(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)+...
            +3/8/f1./u_m_off.*(theta_off(:,3)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_off(:,7)...
            +kappa*theta_off(:,10)).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/f1);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)+...
            +3/8/f1./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_on(:,7)...
            +kappa*theta_on(:,10)).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/f1);
        Phi_e = [Phi_nuoff ; Phi_nuon];
        Phi_pe = gradient(Phi_e)./gradient(u_m_e);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density
        E_psi_e1 = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
        epsilon01 = E_psi_e1(end)*u_m_e(end)^1.5;
        epsilon02 = 3/8*(1-kappa)*kappa21/(1+kappa21)/sqrt(u_meq);
        E_CDM_e = epsilon02.*u_m_e.^-1.5;
        f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);%final interation for beta^2
        delta1 = 1 - E_psi_e1/4;  %delta_omega1 prepared for next stage
        del_psi_e1 = F_e.*theta_e(:,2).*theta_e(:,6)+sin(theta_e(:,1)).*theta_e(:,5)...
            -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,5)./u_m_e;
        Delta_theta_e1 = del_psi_e1./E_psi_e1;
        Delta_theta_e2 = theta_e(:,7);
        Delta_theta_e = (f1*del_psi_e1+E_CDM_e.*Delta_theta_e2)./...
            (f1*E_psi_e1 + E_CDM_e);
        clear del_psi_e1 del_psi_e1
        %%%covariant energy density

        %%%initial condition for Schrodinger
        %%% component 1
        [pk1,loc1] = findpeaks(Delta_theta_e1,u_m_e);
        [pk2,loc2] = findpeaks(-Delta_theta_e1,u_m_e);
        if loc1(end)>loc2(end)
            pk = (pk1(length(loc1)-1)-pk2(end))/2;
        else
            pk = (pk1(length(loc1)-1)-pk2(end-1))/2;
        end
        index_mid = find(Delta_theta_e1==pk1(length(loc1)-1),1,'last');
        diff1 = gradient(pk1)./gradient(loc1);
        diff2 = gradient(pk2)./gradient(loc2);
        %{
            We find that some low k mode would generate incorrect slope since
            findpeaks function would sometimes find both peaks and crests but 
            we want peaks only. Since gradient function subtracts nearby 2 
            points, it's better to use the penultimate point which guarantees 
            it subtracts both peaks or both crests.
        %}
        Delta_theta_e_p1 = 0.5*(diff1(end-1)-diff2(end-1));
        u_mmid = loc1(length(loc1)-1);
        F_emid_sq = sqrt(F_e(index_mid));

        Sc = 3*F_emid_sq/2/u_mmid*(3/32/delta1(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/2/F_emid_sq/u_mmid/delta1(index_mid)-...
            E_psi_e1(index_mid)/4/F_emid_sq/delta1(index_mid));
        Ss = 3/16/delta1(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/F_emid_sq/u_mmid/delta1(index_mid)+3*F_emid_sq/2/u_mmid^2-...
            3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1);
        Dp = Delta_theta_e_p1-3/2/u_mmid/delta1(index_mid)*Phi_e(index_mid)-...
            4*Phi_pe(index_mid)/delta1(index_mid);
        phis01 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid*Sc);
        phic01 = 0.5*(pk+3*F_emid_sq/2/u_mmid*phis01);

        clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p delta1 
        clear Sc Ss Dp
        %%%initial condition for Schrodinger

        %%%reshape exact solution to indexmid
        u_m_e = u_m_e(1:index_mid);
        theta_e = theta_e(1:index_mid,:);
        F_e = F_e(1:index_mid);
        Phi_e = Phi_e(1:index_mid);
        Delta_theta_e1 = Delta_theta_e1(1:index_mid);
        Delta_theta_e2 = Delta_theta_e2(1:index_mid);
        Delta_theta_e = Delta_theta_e(1:index_mid);
        del_ph_p = gradient(theta_e(:,3))./gradient(u_m_e);
        del_b_p = gradient(theta_e(:,10))./gradient(u_m_e);
        %%%reshape exact solution to indexmid


        %%%ode for weak approximation
        [u_m_a,theta_a] = ode45(@(u_m,theta)...
            calculate_cdmtheta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
            alpha_nu,epsilon01,f1,kappa,kappa21),...
            [u_mmid u_mf],...
            [theta_e(end,3);del_ph_p(end);phic01;phis01;...
            theta_e(end,7);theta_e(end,8)...
            ;theta_e(end,9);theta_e(end,10);del_b_p(end)]);
        clear del_ph_p del_b_p
        %%%ode for weak approximation

        %%%covariant energy density
        F_a = 1 + (kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m_a/A_m)+...
            8/3*f1*sqrt(u_m_a)*epsilon01;
        Delta_theta_a1 = 2.*theta_a(:,3)-3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
        Delta_theta_a2 = theta_a(:,5);
        Delta_theta_a = (f1*epsilon01*u_m_a.^-1.5.*Delta_theta_a1+...
            epsilon02*u_m_a.^-1.5.*Delta_theta_a2)./...
            (f1*epsilon01*u_m_a.^-1.5+epsilon02*u_m_a.^-1.5);
        %%%covariant energy density

        %%%Phi for approx
        if theta_e(index_mid,9) == 1 %if it's not turn on before mid
            index_nuturn = find(theta_a(:,7) == 1,1,'last');
            theta_off = theta_a(1:index_nuturn,:);
            u_m_off = u_m_a(1:index_nuturn,:);
            F_off = F_a(1:index_nuturn,:);
            theta_on = theta_a(index_nuturn+1:end,:);
            u_m_on = u_m_a(index_nuturn+1:end,:);
            F_on = F_a(index_nuturn+1:end,:);

            Phi_off = f1*u_m_off.^(-1/2).*epsilon01./A_k.*...
                (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
                3/8/A_k./u_m_off.*(theta_off(:,1)+...
                kappa*theta_off(:,8).*sqrt(u_m_off/A_m)+...
                (1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_off/A_m).*theta_off(:,5));
            Phi_on = f1*u_m_on.^(-1/2)*epsilon01./A_k.*...
                (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
                3/8/A_k./u_m_on.*...
                (theta_on(:,1)./(1+alpha_nu)+...
                kappa*theta_on(:,8).*sqrt(u_m_on/A_m)...
                +(1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_on/A_m).*theta_on(:,5));
            Phi_a = [Phi_off ; Phi_on];
            clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        else
            Phi_a = f1*u_m_a.^(-1/2).*epsilon01/A_k.*...
                (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
                3/8/A_k./u_m_a.*...
                (theta_a(:,1)/(1+alpha_nu)+...
                kappa*theta_a(:,8).*sqrt(u_m_a/A_m)+...
                (1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_a/A_m).*theta_a(:,5));
        end
        %%%Phi for approx

        %%%velocity
        g_theta_a1 = 1.5*sqrt(F_a(end))*theta_a(end,4)/u_m_a(end);
        g_theta_a2 = theta_a(end,6);
        g_theta_a = (f1*epsilon01*g_theta_a1(end)...
            +epsilon02*g_theta_a2(end))./(f1*epsilon01+epsilon02);
        D_r = theta_a(end,2)-theta_a(end,1)/2/u_m_a(end);
        g_r_a = D_r/(2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_m_a(end));
        g_b_a = theta_a(end,9)/(2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_m_a(end));
        %%%velocity
    
        %%%output data
        Delta_theta = Delta_theta_a(end);
        Delta_theta1 = Delta_theta_a1(end);
        Delta_theta2 = Delta_theta_a2(end);
        g_theta = g_theta_a;
        g_theta1 = g_theta_a1;
        g_theta2 = g_theta_a2;
        Delta_b = theta_a(end,8);
        g_b = g_b_a;
        Delta_ph = theta_a(end,1);
        g_ph = g_r_a;
        %%%output data
    else
       %%%theta ode exact
        opts = odeset('MaxStep',res);
        [u_m_e,theta_e] = ode45(@(u_m,theta)...
            calculate_cdmtheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
            f1,kappa,kappa21),[u_m0 u_mmid],...
            [Theta_i1;Theta_ip1;Del_ph0;g_ph0;...
            psi01;psi11;Del_D0;g_D0;nu_turnon],opts);
        Del_b = 3/4*theta_e(:,3);
        g_b = 3/4*theta_e(:,4);
        theta_e = [theta_e Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta_e(:,9) == 1,1,'last');
        F_e = (3/4+3/4*(kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m_e./A_m)+...
            2*f1*u_m_e.^2.*(1-cos(theta_e(:,1))))...
            ./(0.75-f1*u_m_e.^2.*theta_e(:,2).^2);
        theta_off = theta_e(1:index_nuturn,:);
        u_m_off = u_m_e(1:index_nuturn,:);
        F_off = F_e(1:index_nuturn,:);
        theta_on = theta_e(index_nuturn+1:end,:);
        u_m_on = u_m_e(index_nuturn+1:end,:);
        F_on = F_e(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)+...
            +3/8/f1./u_m_off.*(theta_off(:,3)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_off(:,7)...
            +kappa*theta_off(:,10)).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/f1);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)+...
            +3/8/f1./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +((1-kappa)*kappa21/(1+kappa21)*theta_on(:,7)...
            +kappa*theta_on(:,10)).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/f1);
        Phi_e = [Phi_nuoff ; Phi_nuon];
        Phi_pe = gradient(Phi_e)./gradient(u_m_e);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density
        E_psi_e1 = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
        epsilon01 = E_psi_e1(end)*u_m_e(end)^1.5;
        epsilon02 = 3/8*(1-kappa)*kappa21/(1+kappa21)/sqrt(u_meq);
        E_CDM_e = epsilon02.*u_m_e.^-1.5;
        f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);%final interation for beta^2
        delta1 = 1 - E_psi_e1/4;  %delta_omega1 prepared for next stage
        del_psi_e1 = F_e.*theta_e(:,2).*theta_e(:,6)+sin(theta_e(:,1)).*theta_e(:,5)...
            -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,5)./u_m_e;
        Delta_theta_e1 = del_psi_e1./E_psi_e1;
        Delta_theta_e2 = theta_e(:,7);
        Delta_theta_e = (f1*del_psi_e1+E_CDM_e.*Delta_theta_e2)./...
            (f1*E_psi_e1 + E_CDM_e);
        clear del_psi_e1 del_psi_e1
        %%%covariant energy density

        %%%initial condition for Schrodinger
        %%% component 1
        [pk1,loc1] = findpeaks(Delta_theta_e1,u_m_e);
        [pk2,loc2] = findpeaks(-Delta_theta_e1,u_m_e);
        if loc1(end)>loc2(end)
            pk = (pk1(length(loc1)-1)-pk2(end))/2;
        else
            pk = (pk1(length(loc1)-1)-pk2(end-1))/2;
        end
        index_mid = find(Delta_theta_e1==pk1(length(loc1)-1),1,'last');
        diff1 = gradient(pk1)./gradient(loc1);
        diff2 = gradient(pk2)./gradient(loc2);
        %{
            We find that some low k mode would generate incorrect slope since
            findpeaks function would sometimes find both peaks and crests but 
            we want peaks only. Since gradient function subtracts nearby 2 
            points, it's better to use the penultimate point which guarantees 
            it subtracts both peaks or both crests.
        %}
        Delta_theta_e_p1 = 0.5*(diff1(end-1)-diff2(end-1));
        u_mmid = loc1(length(loc1)-1);
        F_emid_sq = sqrt(F_e(index_mid));

        Sc = 3*F_emid_sq/2/u_mmid*(3/32/delta1(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/2/F_emid_sq/u_mmid/delta1(index_mid)-...
            E_psi_e1(index_mid)/4/F_emid_sq/delta1(index_mid));
        Ss = 3/16/delta1(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/F_emid_sq/u_mmid/delta1(index_mid)+3*F_emid_sq/2/u_mmid^2-...
            3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1);
        Dp = Delta_theta_e_p1-3/2/u_mmid/delta1(index_mid)*Phi_e(index_mid)-...
            4*Phi_pe(index_mid)/delta1(index_mid);
        phis01 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid*Sc);
        phic01 = 0.5*(pk+3*F_emid_sq/2/u_mmid*phis01);

        clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p delta1 
        clear Sc Ss Dp
        %%%initial condition for Schrodinger

        %%%reshape exact solution to indexmid
        u_m_e = u_m_e(1:index_mid);
        theta_e = theta_e(1:index_mid,:);
        F_e = F_e(1:index_mid);
        Phi_e = Phi_e(1:index_mid);
        Delta_theta_e1 = Delta_theta_e1(1:index_mid);
        Delta_theta_e2 = Delta_theta_e2(1:index_mid);
        Delta_theta_e = Delta_theta_e(1:index_mid);
        del_ph_p = gradient(theta_e(:,3))./gradient(u_m_e);
        del_b_p = gradient(theta_e(:,10))./gradient(u_m_e);
        %%%reshape exact solution to indexmid


        %%%ode for weak approximation
        [u_m_a,theta_a] = ode45(@(u_m,theta)...
            calculate_cdmtheta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
            alpha_nu,epsilon01,f1,kappa,kappa21),...
            [u_mmid u_mtr],...
            [theta_e(end,3);del_ph_p(end);phic01;phis01;...
            theta_e(end,7);theta_e(end,8)...
            ;theta_e(end,9);theta_e(end,10);del_b_p(end)]);
        clear del_ph_p del_b_p
        %%%ode for weak approximation

        %%%covariant energy density
        F_a = 1 + (kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m_a/A_m)+...
            8/3*f1*sqrt(u_m_a)*epsilon01;
        Delta_theta_a1 = 2.*theta_a(:,3)-3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
        Delta_theta_a2 = theta_a(:,5);
        Delta_theta_a = (f1*epsilon01*u_m_a.^-1.5.*Delta_theta_a1+...
            epsilon02*u_m_a.^-1.5.*Delta_theta_a2)./...
            (f1*epsilon01*u_m_a.^-1.5+epsilon02*u_m_a.^-1.5);
        %%%covariant energy density

        %%%Phi for approx
        if theta_e(index_mid,9) == 1 %if it's not turn on before mid
            index_nuturn = find(theta_a(:,7) == 1,1,'last');
            theta_off = theta_a(1:index_nuturn,:);
            u_m_off = u_m_a(1:index_nuturn,:);
            F_off = F_a(1:index_nuturn,:);
            theta_on = theta_a(index_nuturn+1:end,:);
            u_m_on = u_m_a(index_nuturn+1:end,:);
            F_on = F_a(index_nuturn+1:end,:);

            Phi_off = f1*u_m_off.^(-1/2).*epsilon01./A_k.*...
                (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
                3/8/A_k./u_m_off.*(theta_off(:,1)+...
                kappa*theta_off(:,8).*sqrt(u_m_off/A_m)+...
                (1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_off/A_m).*theta_off(:,5));
            Phi_on = f1*u_m_on.^(-1/2)*epsilon01./A_k.*...
                (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
                3/8/A_k./u_m_on.*...
                (theta_on(:,1)./(1+alpha_nu)+...
                kappa*theta_on(:,8).*sqrt(u_m_on/A_m)...
                +(1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_on/A_m).*theta_on(:,5));
            Phi_a = [Phi_off ; Phi_on];
            clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        else
            Phi_a = f1*u_m_a.^(-1/2).*epsilon01/A_k.*...
                (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
                3/8/A_k./u_m_a.*...
                (theta_a(:,1)/(1+alpha_nu)+...
                kappa*theta_a(:,8).*sqrt(u_m_a/A_m)+...
                (1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_a/A_m).*theta_a(:,5));
        end
        %%%Phi for approx

        %%%reshape to 3rd phase
        [pkr,locr] = findpeaks(abs(theta_a(:,1)),u_m_a);
        index_tr = find(abs(theta_a(:,1))==pkr(length(locr)),1,'last');
        u_m_a = u_m_a(1:index_tr);
        F_a = F_a(1:index_tr);
        theta_a = theta_a(1:index_tr,:);
        Delta_theta_a = Delta_theta_a(1:index_tr);
        Delta_theta_a1 = Delta_theta_a1(1:index_tr);
        Delta_theta_a2 = Delta_theta_a2(1:index_tr);
        Phi_a = Phi_a(1:index_tr);
        u_mtr = u_m_a(index_tr);
        %%%reshape to 3rd phase

        %%%initial condition for 3rd phase
        D_r0 = theta_a(end,2)-theta_a(end,1)/2/u_m_a(end);

        T = T_CMB*sqrt(u_mn/u_mtr);
        alpha_e = 9.78/137^2/m_e^2*sqrt(BE/T)*log(BE/T);
        beta_e = alpha_e*(m_e*T/2/pi)^1.5*exp(-BE/T);
        n_b =  n_b0*(u_mn/u_mtr)^1.5;
        P = beta_e/alpha_e/n_b;
        X = (-P+sqrt(P^2+4*P))/2;
        if abs(X) > 1
            X = 1;
        end
        A_tr = A_t*X;

        Q0 = sqrt(F_a(end))/(4*alpha0/sqrt(u_mtr)+4*beta0/u_mtr+...
            sqrt(F_a(end))/2)*theta_a(end,1);
        g_r0 = (D_r0-A_tr*u_mtr^-1.5*sqrt(F_a(end))*Q0)/...
            (2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_mtr);
        %%%initial condition for 3rd phase

        %%%ode 3rd phase
        [u_m_3,theta_3] = ode45(...
            @(u_m,theta) calculate_cdmtheta_weak_re(u_m,theta,A_k,A_m,A_t,...
            alpha_b,alpha_nu,epsilon01,f1,m_e,BE,T_CMB,u_mn,n_b0,kappa,kappa21),...
            [u_mtr u_mf],...
            [theta_a(end,1);g_r0;theta_a(end,3);theta_a(end,4)...
            ;theta_a(end,5);theta_a(end,6);theta_a(end,7);theta_a(end,8);Q0]);
        %%%ode 3rd phase

        %%%some parameters
        F_3 = 1 + (kappa+(1-kappa)*kappa21/(1+kappa21))*sqrt(u_m_3/A_m)+...
            8./3.*sqrt(u_m_3)*f1*epsilon01;
        Phi_3 = -epsilon01*f1/A_k.*u_m_3.^(-1/2).*(2.*theta_3(:,3)...
            -3.*sqrt(F_3)/2./u_m_3.*theta_3(:,4))...
            -3/8/A_k./u_m_3.*(theta_3(:,1)/(1+alpha_nu)+...
            kappa*sqrt(u_m_3/A_m).*theta_3(:,8)+...
            (1-kappa)*kappa21/(1+kappa21)*sqrt(u_m_3/A_m).*theta_3(:,5));
        %%%some parameters

        %%%covariant energy density 3
        Delta_theta_31 = 2.*theta_3(:,3)-sqrt(F_3)*3/2./u_m_3.*theta_3(:,4);
        Delta_theta_32 = theta_3(:,5);
        Delta_theta_3 = (f1*epsilon01*u_m_3.^-1.5.*Delta_theta_31+...
            epsilon02*u_m_3.^-1.5.*Delta_theta_32)./...
            (f1*epsilon01*u_m_3.^-1.5+epsilon02*u_m_3.^-1.5);
        g_theta_31 = 1.5*sqrt(F_3).*theta_3(:,4)./u_m_3;
        g_theta_32 = theta_3(:,6);
        g_theta_3 = (f1*epsilon01*u_m_3.^-1.5.*g_theta_31...
            +epsilon02*u_m_3.^-1.5.*g_theta_32)./...
            (f1*epsilon01*u_m_3.^-1.5+epsilon02*u_m_3.^-1.5);
        %%%covariant energy density 3

        %%%output data
        Delta_theta = Delta_theta_3(end);
        Delta_theta1 = Delta_theta_31(end);
        Delta_theta2 = Delta_theta_32(end);
        Delta_b = theta_3(end,8);
        Delta_ph = theta_3(end,1);
        g_theta = g_theta_3(end);
        g_theta1 = g_theta_31(end);
        g_theta2 = g_theta_32(end);
        g_b = 0.75*(theta_3(end,2)+theta_3(end,9));
        g_ph = theta_3(end,2);
        %%%output data

    end
    kn(n) = k;
    Delta_thetan(n) = Delta_theta;
    Delta_theta1n(n) = Delta_theta1;
    Delta_theta2n(n) = Delta_theta2;
    Delta_bn(n) = Delta_b(end);
    g_bn(n)= g_b;
    Delta_phn(n) = Delta_ph;
    g_rn(n) = g_ph;
    g_thetan(n) = g_theta;
    g_theta1n(n) = g_theta1;
    g_theta2n(n) = g_theta2;
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

%%% filename
m1 = num2str(mr*1e-22);
af = num2str(sqrt(u_mf/u_meq));
ki = num2str(k_initial);
kf = num2str(k_end);
angle1 = num2str(180 - angle_shift1);
kappa21str = num2str(kappa21);

filename = ['EACDM_two_m1=' m1 '_angle=' angle1 '_1to' ...
    kappa21str '_' af 'a_eq_from_k_' ki 'to' kf '.txt'];
%%% filename
fid = fopen(filename,'wt');
fprintf(fid,'Hubble constant now = %17.24f km/s/Mpc \n',...
    100*Hubble_constant);
fprintf(fid,'Omega_axion = %17.24f \n',Omega_axion);
fprintf(fid,'Omega_baryon = %17.24f \n',Omega_baryon);
fprintf(fid,'Omega_photon = %17.24f \n',Omega_photon);
fprintf(fid,'Omega_neutrino = %17.24f \n',Omega_neutrino);
fprintf(fid,'critical k_c22 = 1/(%3.5f kpc) \n',1/k_c_kpc/sqrt(mr));

fprintf(fid,...
    ['k/k_c \t Delta_theta \t Delta_theta1 \t Delta_theta2 \t Delta_baryon'...
    '\t Delta_photon \t g_A \t g_A1 \t g_A2 \t g_b \t g_r \n']);
fprintf(fid,['%17.24f \t %17.24f\t %17.24f\t %17.24f\t %17.24f\t %17.24f\t'...
    '%17.24f\t %17.24f\t %17.24f\t %17.24f\t %17.24f\t\n'],...
    [kn;Delta_thetan;Delta_theta1n;Delta_theta2n;Delta_bn;...
    Delta_phn;g_thetan;g_theta1n;g_theta2n;g_bn;g_rn]);
fclose(fid);
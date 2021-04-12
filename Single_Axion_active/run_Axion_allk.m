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

%particle mass in 10^-22 eV
mr = 1;

%background axion angle shift from the top of the potential hill:
angle_shift = 90;

%current axion background ratio
Omega_axion = 0.24; 

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
res = 0.1;

%Maximum resolution for iteration of decay constant determination
res_i = 0.01;

%number of iteration for decay constant determination
N_i = 2;

%iteration end point in a_eq
a_stop = 0.05;

%%%input parameters

%%% parameter conversion
G = 6.708830620767732e-57;
a_eq = (Omega_photon+Omega_neutrino)/(Omega_axion+Omega_baryon);
H0 = 2.133119459810834e-33*Hubble_constant;
Hr_eq = H0*sqrt((Omega_photon+Omega_neutrino)/a_eq^2);
u_meq = mr*1e-22*a_eq/2/Hr_eq;
k_c = sqrt(2*mr*1e-22*a_eq*Hr_eq);
k_c_kpc = k_c/1.97326960277e-5/3.24077929e-22;
u_mn = a_eq^-2*u_meq;
rho_cr = 3*H0^2/8/pi/G;
n_b0 = rho_cr/m_n*Omega_baryon;
A_t = n_b0*sigma_T*sqrt(mr*1e-22)/2/sqrt(2)/a_eq^3/...
    (8*pi*G/3*rho_cr*(Omega_baryon+Omega_axion)/a_eq^3)^0.75;
%%% parameter conversion

tic
%%% initial setting
k = k_initial;
u_m0 = 0.0005;          %starting point
u_mf = a_f^2*u_meq;    %end point
n = 1;                  %parameter for save
%%% initial setting

%%%zeroth order parameter
alpha_nu = Omega_neutrino/Omega_photon;
kappa = Omega_baryon/(Omega_baryon+Omega_axion);
alpha_b = (1.0+alpha_nu)*kappa;
alpha_t = (1.0+alpha_nu);
%%%zeroth order parameter

%%% first order parameter
coef_photon = 2/3;
coef_phi = coef_photon*3/2;
alpha0 = A_t/2;
beta0 = 4*sqrt(u_meq)/3/alpha_b*alpha0;
A_m = u_meq;
%%% first order parameter

%%%zeroth order initial condition
Theta0 = pi - angle_shift/180*pi;
Theta_i = Theta0-sin(Theta0)*u_m0^2/5+...
    u_m0^2.5*sin(Theta0)/sqrt(A_m)*11/75;
Theta_ip = -2*sin(Theta0)*u_m0/5+...
    u_m0^2.5*sin(Theta0)/sqrt(A_m)*11/30;
%%%zeroth order initial condition

%%%calculate beta^2 by passive evolution
u_mstop = (a_stop)^2*u_meq;
[u_m_p,theta_p] = ode45(@(u_m,Theta)...
      Psi(u_m,Theta),...
      [u_m0 u_mstop],[Theta_i;Theta_ip]);
epsilon0 = (0.5*theta_p(end,2)^2+1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
beta_2 = 3*(1-kappa)/8/sqrt(A_m)/epsilon0;
    
%%%interation
int = 1; %number of interation
opts = odeset('MaxStep',res_i);
for i = 1:int
       [u_m_p,theta_p] = ode45(@(u_m,Theta)...
            calculate_Theta(u_m,Theta,A_m,beta_2,kappa),...
            [u_m0 u_mstop],[Theta_i;Theta_ip],opts);
       F = (3/4+3/4*(u_m_p./A_m).^(1/2)+2*beta_2*u_m_p.^2.*(1-cos(theta_p(:,1))))./...
            (0.75-beta_2*u_m_p.^2.*theta_p(:,2).^2);
       epsilon0 = (0.5*F(end)*theta_p(end,2)^2+1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
       beta_2 = 3*(1-kappa)/8/sqrt(A_m)/epsilon0;
end
clear u_m_p theta_p u_mstop epsilon0 F int i
%%%calculate beta^2 by passive evolution

while k<=k_end
    %%%k and Transition point for 3rd
    A_k = k^2/mr;
    u_mmid = 30/A_k;
    if u_mmid <= 100
       u_mmid = 100;
    end
    u_mtr = 2*alpha0*0.1/sqrt(A_k);
    %%%k and Transition point for 3rd
    
    %%%initial condition for exact strong
    u_k0 = sqrt(4/3*A_k*u_m0);
    Del_ph0 = coef_photon*(-cos(u_k0)+sin(u_k0)/u_k0...
        -(0.25 +0.5625*alpha_b)*u_k0^2*sqrt(u_m0/A_m)/3.0); 
    g_ph0 = coef_photon*(sin(u_k0)/u_k0+2*cos(u_k0)/u_k0^2-...
        2*sin(u_k0)/u_k0^3+(1 - 0.75*alpha_b)*sqrt(u_m0/A_m)/6.0);
    psi0 = coef_phi*sin(Theta0)*(4*u_m0^2-2*(8/15+alpha_b*21/20)...
        *u_m0^2.5/sqrt(A_m)-176*A_k*u_m0^3/210)/90;
    psi1 = coef_phi*sin(Theta0)*(8*u_m0-5*(8/15+alpha_b*21/20)*...
       u_m0^1.5/sqrt(A_m)-176*A_k*u_m0^2/70)/90;
    nu_turnon = 1;
    %%%initial condition for exact strong

    if u_mf<=u_mmid
        %%%theta ode exact
        opts = odeset('MaxStep',res);
        [u_m,theta] = ode45(@(u_m,theta)...
            calculate_theta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,beta_2,kappa,coef_photon),...
            [u_m0 u_mf],[Theta_i;Theta_ip;Del_ph0;g_ph0;psi0;psi1;nu_turnon],opts);
        Del_b = 3/4*theta(:,3);
        g_b = 3/4*theta(:,4);
        theta = [theta Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta(:,7) == 1,1,'last');
        F = (3/4+3/4*kappa*sqrt(u_m./A_m)+2*beta_2*u_m.^2.*(1-cos(theta(:,1))))...
            ./(0.75-beta_2*u_m.^2.*theta(:,2).^2);
        theta_off = theta(1:index_nuturn,:);
        u_m_off = u_m(1:index_nuturn,:);
        F_off = F(1:index_nuturn,:);
        theta_on = theta(index_nuturn+1:end,:);
        u_m_on = u_m(index_nuturn+1:end,:);
        F_on = F(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)...
            +3/8/beta_2./u_m_off.*(theta_off(:,3)...
            +kappa*theta_off(:,8).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/beta_2);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)...
            +3/8/beta_2./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +kappa*theta_on(:,8).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/beta_2);
        Phi = [Phi_nuoff ; Phi_nuon];
        Phi_p = gradient(Phi)./gradient(u_m);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density
        E_psi = 0.5*F.*theta(:,2).^2+1-cos(theta(:,1));
        del_psi = F.*theta(:,2).*theta(:,6)+sin(theta(:,1)).*theta(:,5)-Phi.*F.*theta(:,2).^2+...
            1.5.*F.*theta(:,2).*theta(:,5)./u_m;
        Delta_theta = del_psi./E_psi;
        %%%covariant energy density
        
        %%%final value
        Delta_theta = Delta_theta(end);
        Delta_b = theta(end,8);
        Delta_ph = theta(end,3);
        g_A = 1.5.*F(end).*theta(end,2).*theta(end,5)./u_m(end)/E_psi(end);
        gr = theta(end,4);
        gb = theta(end,9);
        %%%final value

    elseif u_mf<=u_mtr
        %%%theta ode exact
        opts = odeset('MaxStep',res);
        [u_m_e,theta_e] = ode45(@(u_m,theta)...
            calculate_theta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,beta_2,kappa,coef_photon),...
            [u_m0 u_mmid],[Theta_i;Theta_ip;Del_ph0;g_ph0;psi0;psi1;nu_turnon],opts);
        Del_b = 3/4*theta_e(:,3);
        g_b = 3/4*theta_e(:,4);
        theta_e = [theta_e Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta_e(:,7) == 1,1,'last');
        F_e = (3/4+3/4*kappa*sqrt(u_m_e./A_m)+2*beta_2*u_m_e.^2.*(1-cos(theta_e(:,1))))...
            ./(0.75-beta_2*u_m_e.^2.*theta_e(:,2).^2);
        theta_off = theta_e(1:index_nuturn,:);
        u_m_off = u_m_e(1:index_nuturn,:);
        F_off = F_e(1:index_nuturn,:);
        theta_on = theta_e(index_nuturn+1:end,:);
        u_m_on = u_m_e(index_nuturn+1:end,:);
        F_on = F_e(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)...
            +3/8/beta_2./u_m_off.*(theta_off(:,3)...
            +kappa*theta_off(:,8).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/beta_2);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)...
            +3/8/beta_2./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +kappa*theta_on(:,8).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/beta_2);
        Phi_e = [Phi_nuoff ; Phi_nuon];
        Phi_pe = gradient(Phi_e)./gradient(u_m_e);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density
        E_psi_e = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
        epsilon0 = E_psi_e(end)*u_m_e(end)^1.5;
        beta_2 = 3*(1-kappa)/8/sqrt(A_m)/epsilon0; %2nd interation for beta^2
        delta = 1 - E_psi_e/4;  %delta_omega prepared for next stage
        del_psi_e = F_e.*theta_e(:,2).*theta_e(:,6)+sin(theta_e(:,1)).*theta_e(:,5)...
            -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,5)./u_m_e;
        Delta_theta_e = del_psi_e./E_psi_e;
        clear del_psi_e 
        %%%covariant energy density

        %%%initial condition for Schrodinger
        [pk1,loc1] = findpeaks(Delta_theta_e,u_m_e);
        [pk2,loc2] = findpeaks(-Delta_theta_e,u_m_e);
        if loc1(end)>loc2(end)
            pk = (pk1(length(loc1)-1)-pk2(end))/2;
        else
            pk = (pk1(length(loc1)-1)-pk2(end-1))/2;
        end
        index_mid = find(Delta_theta_e==pk1(length(loc1)-1),1,'last');
        diff1 = gradient(pk1)./gradient(loc1);
        diff2 = gradient(pk2)./gradient(loc2);
        %{
            We find that some low k mode would generate incorrect slope since
            findpeaks function would sometimes find both peaks and crests but 
            we want peaks only. Since gradient function subtracts nearby 2 
            points, it's better to use the penultimate point which guarantees 
            it subtracts both peaks or both crests.
        %}
        Delta_theta_e_p = 0.5*(diff1(end-1)-diff2(end-1));
        if k^2<=0.004
            Delta_theta_e_p = (Delta_theta_e(end)-Delta_theta_e(end-1000))/...
                (u_m_e(end)-u_m_e(end-1000));
        end
        u_mmid = loc1(length(loc1)-1);
        F_emid_sq = sqrt(F_e(index_mid));

        Sc = 3*F_emid_sq/2/u_mmid*(3/32/delta(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/2/F_emid_sq/u_mmid/delta(index_mid)-...
            E_psi_e(index_mid)/4/F_emid_sq/delta(index_mid));
        Ss = 3/16/delta(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/F_emid_sq/u_mmid/delta(index_mid)+3*F_emid_sq/2/u_mmid^2-...
            3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1);
        Dp = Delta_theta_e_p-3/2/u_mmid/delta(index_mid)*Phi_e(index_mid)-...
            4*Phi_pe(index_mid)/delta(index_mid);
        phis0 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid*Sc);
        phic0 = 0.5*(pk+3*F_emid_sq/2/u_mmid*phis0);

        clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p F_emid_sq ...
            delta Phi_pe
        clear Sc Ss Dp
        %%%initial condition for Schrodinger

        %%%reshape exact solution to indexmid
        u_m_e = u_m_e(1:index_mid);
        theta_e = theta_e(1:index_mid,:);
        F_e = F_e(1:index_mid);
        Phi_e = Phi_e(1:index_mid);
        Delta_theta_e = Delta_theta_e(1:index_mid);
        del_ph_p = gradient(theta_e(:,3))./gradient(u_m_e);
        del_b_p = gradient(theta_e(:,8))./gradient(u_m_e);
        %%%reshape exact solution to indexmid


        %%%ode for weak approximation
        [u_m_a,theta_a] = ode45(@(u_m,theta)...
            calculate_theta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
        alpha_nu,epsilon0,beta_2,coef_photon,kappa),...
            [u_mmid u_mf],...
            [theta_e(end,3);del_ph_p(end);phic0;phis0...
            ;theta_e(end,7);theta_e(end,8);del_b_p(end)]);
        clear del_ph_p del_b_p
        %%%ode for weak approximation

        %%%covariant energy density
        F_a = 1 + kappa*sqrt(u_m_a/A_m)+8./3.*beta_2*sqrt(u_m_a)*epsilon0;
        Delta_theta_a = 2.*theta_a(:,3)...
            -3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
        %%%covariant energy density

        %%%Phi for approx
        if theta_e(index_mid,7) == 1 %if it's not turn on before mid
            index_nuturn = find(theta_a(:,5) == 1,1,'last');
            theta_off = theta_a(1:index_nuturn,:);
            u_m_off = u_m_a(1:index_nuturn,:);
            F_off = F_a(1:index_nuturn,:);
            theta_on = theta_a(index_nuturn+1:end,:);
            u_m_on = u_m_a(index_nuturn+1:end,:);
            F_on = F_a(index_nuturn+1:end,:);

            Phi_off = beta_2*u_m_off.^(-1/2).*epsilon0./A_k.*...
                (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
                3/8/A_k./u_m_off.*(theta_off(:,1)+kappa*theta_off(:,6).*sqrt(u_m_off/A_m));
            Phi_on = beta_2*u_m_on.^(-1/2)*epsilon0./A_k.*...
                (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
                3/8/A_k./u_m_on.*...
                (theta_on(:,1)./(1+alpha_nu)+kappa*theta_on(:,6).*sqrt(u_m_on/A_m));
            Phi_a = [Phi_off ; Phi_on];
            clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        else
            Phi_a = beta_2*u_m_a.^(-1/2).*epsilon0/A_k.*...
                (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
                3/8/A_k./u_m_a.*...
                (theta_a(:,1)/(1+alpha_nu)+kappa*theta_a(:,6).*sqrt(u_m_a/A_m));
        end
        %%%Phi for approx

        %%%reshape all
        Delta_theta = Delta_theta_a(end);
        Delta_b = theta_a(end,6);
        Delta_ph = theta_a(end,1);
        g_A = 1.5*sqrt(F_a(end))*theta_a(end,4)/u_m_a(end);
        D_r0 = theta_a(end,2)-theta_a(end,1)/2/u_m_a(end);
        gr = (D_r0)/(2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_m_a(end));
        gb = theta_a(end,7)/(2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_m_a(end));
        %%%reshape all
    else
        %%%theta ode exact
        opts = odeset('MaxStep',0.1);
        [u_m_e,theta_e] = ode45(@(u_m,theta)...
            calculate_theta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,beta_2,kappa,coef_photon),...
            [u_m0 u_mmid],[Theta_i;Theta_ip;Del_ph0;g_ph0;psi0;psi1;nu_turnon],opts);
        Del_b = 3/4*theta_e(:,3);
        g_b = 3/4*theta_e(:,4);
        theta_e = [theta_e Del_b g_b];
        clear Del_b g_b
        %%%theta ode exact

        %%% calculate phi
        index_nuturn = find(theta_e(:,7) == 1,1,'last');
        F_e = (3/4+3/4*kappa*sqrt(u_m_e./A_m)+2*beta_2*u_m_e.^2.*(1-cos(theta_e(:,1))))...
            ./(0.75-beta_2*u_m_e.^2.*theta_e(:,2).^2);
        theta_off = theta_e(1:index_nuturn,:);
        u_m_off = u_m_e(1:index_nuturn,:);
        F_off = F_e(1:index_nuturn,:);
        theta_on = theta_e(index_nuturn+1:end,:);
        u_m_on = u_m_e(index_nuturn+1:end,:);
        F_on = F_e(index_nuturn+1:end,:);

        Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,6)+...
            u_m_off.*sin(theta_off(:,1)).*theta_off(:,5)+...
            1.5.*F_off.*theta_off(:,2).*theta_off(:,5)...
            +3/8/beta_2./u_m_off.*(theta_off(:,3)...
            +kappa*theta_off(:,8).*sqrt(u_m_off/A_m)))./...
            (u_m_off.*F_off.*theta_off(:,2).^2-A_k/beta_2);
        Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,6)...
            +u_m_on.*sin(theta_on(:,1)).*theta_on(:,5)+...
            1.5.*F_on.*theta_on(:,2).*theta_on(:,5)...
            +3/8/beta_2./u_m_on.*(theta_on(:,3)/(1+alpha_nu)...
            +kappa*theta_on(:,8).*sqrt(u_m_on/A_m)))...
            ./(u_m_on.*F_on.*theta_on(:,2).^2-A_k/beta_2);
        Phi_e = [Phi_nuoff ; Phi_nuon];
        Phi_pe = gradient(Phi_e)./gradient(u_m_e);
        clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        %%% calculate phi

        %%%covariant energy density
        E_psi_e = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
        epsilon0 = E_psi_e(end)*u_m_e(end)^1.5;
        beta_2 = 3*(1-kappa)/8/sqrt(A_m)/epsilon0; %2nd interation for beta^2
        delta = 1 - E_psi_e/4;  %delta_omega prepared for next stage
        del_psi_e = F_e.*theta_e(:,2).*theta_e(:,6)+sin(theta_e(:,1)).*theta_e(:,5)...
            -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,5)./u_m_e;
        Delta_theta_e = del_psi_e./E_psi_e;
        clear del_psi_e
        %%%covariant energy density

        %%%initial condition for Schrodinger
        [pk1,loc1] = findpeaks(Delta_theta_e,u_m_e);
        [pk2,loc2] = findpeaks(-Delta_theta_e,u_m_e);
        if loc1(end)>loc2(end)
            pk = (pk1(length(loc1)-1)-pk2(end))/2;
        else
            pk = (pk1(length(loc1)-1)-pk2(end-1))/2;
        end
        index_mid = find(Delta_theta_e==pk1(length(loc1)-1),1,'last');
        diff1 = gradient(pk1)./gradient(loc1);
        diff2 = gradient(pk2)./gradient(loc2);
        %{
            We find that some low k mode would generate incorrect slope since
            findpeaks function would sometimes find both peaks and crests but 
            we want peaks only. Since gradient function subtracts nearby 2 
            points, it's better to use the penultimate point which guarantees 
            it subtracts both peaks or both crests.
        %}
        Delta_theta_e_p = 0.5*(diff1(end-1)-diff2(end-1));
        if k^2<=0.004
            Delta_theta_e_p = (Delta_theta_e(end)-Delta_theta_e(end-1000))/...
                (u_m_e(end)-u_m_e(end-1000));
        end
        u_mmid = loc1(length(loc1)-1);
        F_emid_sq = sqrt(F_e(index_mid));

        Sc = 3*F_emid_sq/2/u_mmid*(3/32/delta(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/2/F_emid_sq/u_mmid/delta(index_mid)-...
            E_psi_e(index_mid)/4/F_emid_sq/delta(index_mid));
        Ss = 3/16/delta(index_mid)/F_emid_sq/u_mmid^2+...
            A_k/F_emid_sq/u_mmid/delta(index_mid)+3*F_emid_sq/2/u_mmid^2-...
            3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1);
        Dp = Delta_theta_e_p-3/2/u_mmid/delta(index_mid)*Phi_e(index_mid)-...
            4*Phi_pe(index_mid)/delta(index_mid);
        phis0 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid*Sc);
        phic0 = 0.5*(pk+3*F_emid_sq/2/u_mmid*phis0);


        clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p F_emid_sq ...
            delta Phi_pe
        clear Sc Ss Dp

        %%%initial condition for Schrodinger

        %%%reshape exact solution to indexmid
        u_m_e = u_m_e(1:index_mid);
        theta_e = theta_e(1:index_mid,:);
        F_e = F_e(1:index_mid);
        Phi_e = Phi_e(1:index_mid);
        Delta_theta_e = Delta_theta_e(1:index_mid);
        del_ph_p = gradient(theta_e(:,3))./gradient(u_m_e);
        del_b_p = gradient(theta_e(:,8))./gradient(u_m_e);
        %%%reshape exact solution to indexmid


        %%%ode for weak approximation
        [u_m_a,theta_a] = ode45(@(u_m,theta)...
            calculate_theta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
        alpha_nu,epsilon0,beta_2,coef_photon,kappa),...
            [u_mmid u_mtr],...
            [theta_e(end,3);del_ph_p(end);phic0;phis0...
            ;theta_e(end,7);theta_e(end,8);del_b_p(end)]);
        clear del_ph_p del_b_p
        %%%ode for weak approximation

        %%%covariant energy density
        F_a = 1 + kappa*sqrt(u_m_a/A_m)+8./3.*beta_2*sqrt(u_m_a)*epsilon0;
        Delta_theta_a = 2.*theta_a(:,3)...
            -3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
        %%%covariant energy density

        %%%Phi for approx
        if theta_e(index_mid,7) == 1 %if it's not turn on before mid
            index_nuturn = find(theta_a(:,5) == 1,1,'last');
            theta_off = theta_a(1:index_nuturn,:);
            u_m_off = u_m_a(1:index_nuturn,:);
            F_off = F_a(1:index_nuturn,:);
            theta_on = theta_a(index_nuturn+1:end,:);
            u_m_on = u_m_a(index_nuturn+1:end,:);
            F_on = F_a(index_nuturn+1:end,:);

            Phi_off = beta_2*u_m_off.^(-1/2).*epsilon0./A_k.*...
                (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
                3/8/A_k./u_m_off.*(theta_off(:,1)+kappa*theta_off(:,6).*sqrt(u_m_off/A_m));
            Phi_on = beta_2*u_m_on.^(-1/2)*epsilon0./A_k.*...
                (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
                3/8/A_k./u_m_on.*...
                (theta_on(:,1)./(1+alpha_nu)+kappa*theta_on(:,6).*sqrt(u_m_on/A_m));
            Phi_a = [Phi_off ; Phi_on];
            clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
        else
            Phi_a = beta_2*u_m_a.^(-1/2).*epsilon0/A_k.*...
                (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
                3/8/A_k./u_m_a.*...
                (theta_a(:,1)/(1+alpha_nu)+kappa*theta_a(:,6).*sqrt(u_m_a/A_m));
        end
        %%%Phi for approx

        %%%reshape to 3rd phase
        [pkr,locr] = findpeaks(abs(theta_a(:,1)),u_m_a);
        index_tr = find(abs(theta_a(:,1))==pkr(length(locr)),1,'last');
        u_m_a = u_m_a(1:index_tr);
        F_a = F_a(1:index_tr);
        theta_a = theta_a(1:index_tr,:);
        Delta_theta_a = Delta_theta_a(1:index_tr);
        Phi_a = Phi_a(1:index_tr);
        u_mtr = u_m_a(index_tr);
        %%%reshape to 3rd phase

        %%%initial condition for 3rd phase
        D_r0 = theta_a(end,2)-theta_a(end,1)/2/u_m_a(end);

        T = T_CMB*sqrt(u_mn/u_mtr);
        n_b =  n_b0*(u_mn/u_mtr)^1.5;
        P = (m_e*T/2/pi)^1.5*exp(-BE/T)/n_b;
        X = (-P+sqrt(P^2+4*P))/2;
        if abs(X) > 1
            X = 1;
        end
        A_tr = A_t*X;

        Q0 = sqrt(F_a(end))/(4*alpha0/sqrt(u_mtr)+4*beta0/u_mtr+...
            sqrt(F_a(end))/2)*theta_a(end,1);
        g_r0 = (D_r0-A_tr/u_mtr^1.5/sqrt(F_a(end))*Q0)/...
            (2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_mtr);
        %%%initial condition for 3rd phase

        %%%ode 3rd phase
        [u_m_3,theta_3] = ode45(...
        @(u_m,theta) calculate_theta_weak_re(u_m,theta,A_k,A_m,A_t,alpha_b...
        ,alpha_nu,epsilon0,beta_2,m_e,BE,T_CMB,u_mn,n_b0,kappa),...
            [u_mtr u_mf],...
            [theta_a(end,1);g_r0;theta_a(end,3);theta_a(end,4)...
            ;theta_a(end,5);theta_a(end,6);Q0]);
        %%%ode 3rd phase

        %%%some parameters
        F_3 = 1 + kappa*sqrt(u_m_3/A_m)+8./3.*sqrt(u_m_3)*beta_2*epsilon0;
        Phi_3 = -epsilon0*beta_2/A_k.*u_m_3.^(-1/2).*(2.*theta_3(:,3)...
            -3.*sqrt(F_3)/2./u_m_3.*theta_3(:,4))...
            -3/8/A_k./u_m_3.*(theta_3(:,1)/(1+alpha_nu)+kappa*sqrt(u_m_3/A_m).*theta_3(:,6));
        %%%some parameters

        %%%covariant energy density 3
        Delta_theta_3 = 2.*theta_3(:,3)-sqrt(F_3)*3/2./u_m_3.*theta_3(:,4);
        %%%covariant energy density 3

        %%%reshape all
        Delta_theta = Delta_theta_3(end);
        Delta_b = theta_3(end,6);
        Delta_ph = theta_3(end,1);
        g_A = 1.5*sqrt(F_3(end))*theta_3(end,4)/u_m_3(end);
        gb = 0.75*(theta_3(end,2)+theta_3(end,7));
        gr = theta_3(end,2);
        %%%reshape all

    end
    kn(n) = k;
    Delta_thetan(n) = Delta_theta(end);
    Delta_bn(n) = Delta_b(end);
    g_bn(n)= gb;
    Delta_phn(n) = Delta_ph(end);
    g_rn(n) = gr;
    theta_An(n) = g_A;
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

%%%save file
m = num2str(mr*1e-22);
angle = num2str(180 - angle_shift);
af = num2str(sqrt(u_mf/u_meq));
ki = num2str(k_initial);
kf = num2str(k_end);
filename = ['Axion_' m '_' angle '_' af 'aeq' '_from_k_' ki 'to' kf '.txt'];
fid = fopen(filename,'wt');
fprintf(fid,'Hubble constant now = %17.24f km/s/Mpc \n',...
    100*Hubble_constant);
fprintf(fid,'Omega_axion = %17.24f \n',Omega_axion);
fprintf(fid,'Omega_baryon = %17.24f \n',Omega_baryon);
fprintf(fid,'Omega_photon = %17.24f \n',Omega_photon);
fprintf(fid,'Omega_neutrino = %17.24f \n',Omega_neutrino);
fprintf(fid,'critical k_c = 1/(%3.5f kpc) \n',1/k_c_kpc);
fprintf(fid,'critical k_c22 = 1/(%3.5f kpc) \n',1/k_c_kpc/sqrt(mr));

fprintf(fid,...
    'k/k_c \t Delta_theta \t Delta_baryon \t Delta_photon \t g_A \t g_b \t g_r \n');
fprintf(fid,'%17.24f \t %17.24f\t %17.24f\t %17.24f\t %17.24f\t %17.24f\t %17.24f\t\n',...
    [kn;Delta_thetan;Delta_bn;Delta_phn;theta_An;g_bn;g_rn]);
fclose(fid);
%%%save file
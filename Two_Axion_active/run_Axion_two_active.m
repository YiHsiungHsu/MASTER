clear all
close all

%%%input parameters
%End time in a_eq (scale factor at matter-radiation equality)
a_f = 3;

%Wave mode in critical k_c22 for 10^-22eV (typically k_c22~(100 kpc)^-1)
k = 1;

%first particle mass in 10^-22 eV
mr = 1;

%second particle mass to first particle ratio (>=1)
mu = 10;

%background axion angle shift from the top of the potential hill for particle 1:
angle_shift1 = 90;

%background axion angle shift from the top of the potential hill for particle 2:
angle_shift2 = 0.04;

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
tic
%%% initial setting
A_k = k^2/mr;             %(k/kcritial)^2
u_m0 = 0.0005*mr;          %starting point
u_mf = a_f^2*u_meq;    %end point

%%% One should adjust resolution of time if changing masses once unexpected
%%% result is obtained. I have put a general adjustment but it may not fit
%%% all cases.

%%% One should also note that there may be cases where findpeaks function at
%%% transition from phase1 to 2 doesn't works well. My solution is using a
%%% overall funtion to enhance prominences of peaks. This isn't a permanent
%%% solution. One can change it if one has better solution. See where I set
%%% p at section "Thomson scattering and others".

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
%%% first order parameter

%%%Thomson scattering and others
%Thomson
alpha0 = A_t/2;
beta0 = 4*sqrt(u_meq)/3/alpha_b*alpha0;
%ohters
A_m = u_meq;
u_mmid = 30/A_k;
if u_mmid <= 100
    u_mmid = 100;
end
u_mtr = 2*alpha0*0.1/sqrt(A_k);
% set power so that find peak could work properly at transition
if A_k/mu >=0.9
        p = 1.5;
    elseif A_k/mu <=0.9 && A_k/mu >=0.8
        p = 1.3;
    elseif A_k/mu <=0.01
        p = -1;
    else
        p = 0;
end
%%%Thomeson Scattering and others

%%%initial condition for exact strong
u_k0 = sqrt(4/3*A_k*u_m0);
Theta01 = pi - angle_shift1/180*pi;
Theta02 = pi - angle_shift2/180*pi;

Theta_i1 = Theta01-sin(Theta01)*u_m0^2/5+...
    u_m0^2.5*sin(Theta01)/sqrt(A_m)*11/75;
Theta_ip1 = -2*sin(Theta01)*u_m0/5+...
    u_m0^1.5*sin(Theta01)/sqrt(A_m)*11/30;
Theta_i2 = Theta02-sin(Theta02)*u_m0^2/5*mu^2+...
    u_m0^2.5*sin(Theta02)/sqrt(A_m)*11/75*mu^2;
Theta_ip2 = -2*sin(Theta02)*u_m0/5*mu^2+...
    u_m0^1.5*sin(Theta02)/sqrt(A_m)*11/30*mu^2;
Del_ph0 = coef_photon*(-cos(u_k0)+sin(u_k0)/u_k0...
    -(0.25 +0.5625*alpha_b)*u_k0^2*sqrt(u_m0/A_m)/3.0); 
g_ph0 = coef_photon*(sin(u_k0)/u_k0+2*cos(u_k0)/u_k0^2-...
    2*sin(u_k0)/u_k0^3+(1 - 0.75*alpha_b)*sqrt(u_m0/A_m)/6.0);
psi01 = coef_phi*sin(Theta01)*(4*u_m0^2-2*(8/15+alpha_b*21/20)...
    *u_m0^2.5/sqrt(A_m)-176*A_k*u_m0^3/210)/90;
psi11 = coef_phi*sin(Theta01)*(8*u_m0-5*(8/15+alpha_b*21/20)*...
    u_m0^1.5/sqrt(A_m)-176*A_k*u_m0^2/70)/90;
psi02 = coef_phi*sin(Theta02)*(4*u_m0^2-2*(8/15+alpha_b*21/20)...
    *u_m0^2.5/sqrt(A_m)-176*A_k*u_m0^3/210)/90*mu^2;
psi12 = coef_phi*sin(Theta02)*(8*u_m0-5*(8/15+alpha_b*21/20)*...
    u_m0^1.5/sqrt(A_m)-176*A_k*u_m0^2/70)/90*mu^2;
nu_turnon = 1;
%%%initial condition for exact strong

%%%calculate beta^2 by passive evolution abd Newton method
u_mstop = a_stop^2*u_meq;
[u_m_p,theta_p] = ode45(@(u_m,Theta)...
        calculate_ByTheta_passive(u_m,Theta,mu),...
        [u_m0 u_mstop],[Theta_i1;Theta_ip1;Theta_i2;Theta_ip2]);
epsilon01 = (0.5*theta_p(end,2)^2+1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
epsilon02 = (0.5*theta_p(end,4)^2+mu^2*(1-cos(theta_p(end,3))))*u_m_p(end)^1.5;
f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
f2 = kappa21*f1*epsilon01/epsilon02;
%%%interation
opts = odeset('MaxStep',res_i/mu);
for i = 1:N_i
    [u_m_p,theta_p] = ode45(@(u_m,Theta)...
            calculate_ByTheta(u_m,Theta,A_m,f1,f2,mu,kappa),...
            [u_m0 u_mstop],[Theta_i1;Theta_ip1;Theta_i2;Theta_ip2],opts);
    F = (3/4+3/4*kappa*(u_m_p./A_m).^(1/2)+2*f1*u_m_p.^2.*(1-cos(theta_p(:,1)))...
        +2*mu^2*f2*u_m_p.^2.*(1-cos(theta_p(:,3))))./...
        (0.75-f1*u_m_p.^2.*theta_p(:,2).^2-f2*u_m_p.^2.*theta_p(:,4).^2);
    epsilon01 = (0.5*F(end)*theta_p(end,2)^2+...
        1-cos(theta_p(end,1)))*u_m_p(end)^1.5;
    epsilon02 = (0.5*F(end)*theta_p(end,4)^2+...
        mu^2*(1-cos(theta_p(end,3))))*u_m_p(end)^1.5;
    f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
    f2 = kappa21*f1*epsilon01/epsilon02;
end

clear u_m_p theta_p u_mstop epsilon0 F int i
%%%calculate beta^2 by passive evolution and Newton method

if u_mf<=u_mmid
    %%%theta ode exact
    opts = odeset('MaxStep',res/mu);
    [u_m,theta] = ode45(@(u_m,theta)...
        calculate_bytheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
        f1,f2,mu,kappa,coef_photon),[u_m0 u_mf],...
        [Theta_i1;Theta_ip1;Theta_i2;Theta_ip2;Del_ph0;g_ph0;...
        psi01;psi11;psi02;psi12;nu_turnon],opts);
    Del_b = 3/4*theta(:,5);
    g_b = 3/4*theta(:,6);
    theta = [theta Del_b g_b];
    clear Del_b g_b
    %%%theta ode exact
    
    %%% calculate phi
    index_nuturn = find(theta(:,11) == 1,1,'last');
    F = (3/4+3/4*kappa*sqrt(u_m./A_m)+2*f1*u_m.^2.*(1-cos(theta(:,1)))+...
        2*mu^2*f2*u_m.^2.*(1-cos(theta(:,3))))...
        ./(0.75-f1*u_m.^2.*theta(:,2).^2-f2*u_m.^2.*theta(:,4).^2);
    theta_off = theta(1:index_nuturn,:);
    u_m_off = u_m(1:index_nuturn,:);
    F_off = F(1:index_nuturn,:);
    theta_on = theta(index_nuturn+1:end,:);
    u_m_on = u_m(index_nuturn+1:end,:);
    F_on = F(index_nuturn+1:end,:);
    
    Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,8)+...
        u_m_off.*sin(theta_off(:,1)).*theta_off(:,7)+...
        1.5.*F_off.*theta_off(:,2).*theta_off(:,7)+...
        f2/f1*(u_m_off.*F_off.*theta_off(:,4).*theta_off(:,10)+...
        mu^2*u_m_off.*sin(theta_off(:,3)).*theta_off(:,9)+...
        1.5.*F_off.*theta_off(:,4).*theta_off(:,9))+...
        +3/8/f1./u_m_off.*(theta_off(:,5)...
        +kappa*theta_off(:,12).*sqrt(u_m_off/A_m)))./...
        (u_m_off.*F_off.*(theta_off(:,2).^2+f2/f1*theta_off(:,4).^2)...
        -A_k/f1);
    Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,8)...
        +u_m_on.*sin(theta_on(:,1)).*theta_on(:,7)+...
        1.5.*F_on.*theta_on(:,2).*theta_on(:,7)+...
        f2/f1*(u_m_on.*F_on.*theta_on(:,4).*theta_on(:,10)...
        +mu^2*u_m_on.*sin(theta_on(:,3)).*theta_on(:,9)+...
        1.5.*F_on.*theta_on(:,4).*theta_on(:,9))...
        +3/8/f1./u_m_on.*(theta_on(:,5)/(1+alpha_nu)...
        +kappa*theta_on(:,12).*sqrt(u_m_on/A_m)))...
        ./(u_m_on.*F_on.*(theta_on(:,2).^2+f2/f1*theta_on(:,4).^2)...
        -A_k/f1);
    Phi = [Phi_nuoff ; Phi_nuon];
    Phi_p = gradient(Phi)./gradient(u_m);
    clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
    %%% calculate phi
    
    %%%covariant energy density
    E_psi1 = f1*(0.5*F.*theta(:,2).^2+1-cos(theta(:,1)));
    E_psi2 = f2*(0.5*F.*theta(:,4).^2+mu^2*(1-cos(theta(:,2))));
    E_psi = E_psi1+E_psi2;
    del_psi1 = f1*(F.*theta(:,2).*theta(:,8)+sin(theta(:,1)).*theta(:,7)-...
        Phi.*F.*theta(:,2).^2+1.5.*F.*theta(:,2).*theta(:,7)./u_m);
    del_psi2 = f2*(F.*theta(:,4).*theta(:,10)+mu^2*sin(theta(:,3)).*theta(:,9)-...
        Phi.*F.*theta(:,4).^2+1.5.*F.*theta(:,4).*theta(:,9)./u_m);
    del_psi = del_psi1+del_psi2;
    Delta_theta = del_psi./E_psi;
    Delta_theta1 = del_psi1./E_psi1;
    Delta_theta2 = del_psi2./E_psi2;
    Delta_b = theta(:,12);
    Delta_ph = theta(:,5);
    g_theta1 = f1*1.5.*F.*theta(:,2).*theta(:,7)./u_m./E_psi1;
    g_theta2 = f2*1.5.*F.*theta(:,4).*theta(:,9)./u_m./E_psi2;
    g_theta = (E_psi1.*g_theta1+E_psi2.*g_theta2)./E_psi;
    g_b = theta(:,13);
    g_ph = theta(:,6);
    %%%covariant energy density
elseif u_mf<=u_mtr
    %%%theta ode exact
    opts = odeset('MaxStep',res/mu);
    [u_m_e,theta_e] = ode45(@(u_m,theta)...
        calculate_bytheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
        f1,f2,mu,kappa,coef_photon),[u_m0 u_mmid],...
        [Theta_i1;Theta_ip1;Theta_i2;Theta_ip2;Del_ph0;g_ph0;...
        psi01;psi11;psi02;psi12;nu_turnon],opts);
    Del_b = 3/4*theta_e(:,5);
    g_b = 3/4*theta_e(:,6);
    theta_e = [theta_e Del_b g_b];
    clear Del_b g_b
    %%%theta ode exact
    
    %%% calculate phi
    index_nuturn = find(theta_e(:,11) == 1,1,'last');
    F_e = (3/4+3/4*kappa*sqrt(u_m_e./A_m)+2*f1*u_m_e.^2.*(1-cos(theta_e(:,1)))...
        +mu^2*2*f2*u_m_e.^2.*(1-cos(theta_e(:,3))))...
        ./(0.75-f1*u_m_e.^2.*theta_e(:,2).^2-f2*u_m_e.^2.*theta_e(:,4).^2);
    theta_off = theta_e(1:index_nuturn,:);
    u_m_off = u_m_e(1:index_nuturn,:);
    F_off = F_e(1:index_nuturn,:);
    theta_on = theta_e(index_nuturn+1:end,:);
    u_m_on = u_m_e(index_nuturn+1:end,:);
    F_on = F_e(index_nuturn+1:end,:);
    
    Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,8)+...
        u_m_off.*sin(theta_off(:,1)).*theta_off(:,7)+...
        1.5.*F_off.*theta_off(:,2).*theta_off(:,7)+...
        f2/f1*(u_m_off.*F_off.*theta_off(:,4).*theta_off(:,10)+...
        mu^2*u_m_off.*sin(theta_off(:,3)).*theta_off(:,9)+...
        1.5.*F_off.*theta_off(:,4).*theta_off(:,9))+...
        +3/8/f1./u_m_off.*(theta_off(:,5)...
        +kappa*theta_off(:,12).*sqrt(u_m_off/A_m)))./...
        (u_m_off.*F_off.*(theta_off(:,2).^2+f2/f1*theta_off(:,4).^2)...
        -A_k/f1);
    Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,8)...
        +u_m_on.*sin(theta_on(:,1)).*theta_on(:,7)+...
        1.5.*F_on.*theta_on(:,2).*theta_on(:,7)+...
        f2/f1*(u_m_on.*F_on.*theta_on(:,4).*theta_on(:,10)...
        +mu^2*u_m_on.*sin(theta_on(:,3)).*theta_on(:,9)+...
        1.5.*F_on.*theta_on(:,4).*theta_on(:,9))...
        +3/8/f1./u_m_on.*(theta_on(:,5)/(1+alpha_nu)...
        +kappa*theta_on(:,12).*sqrt(u_m_on/A_m)))...
        ./(u_m_on.*F_on.*(theta_on(:,2).^2+f2/f1*theta_on(:,4).^2)...
        -A_k/f1);
    Phi_e = [Phi_nuoff ; Phi_nuon];
    Phi_pe = gradient(Phi_e)./gradient(u_m_e);
    clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
    %%% calculate phi
    
    %%%covariant energy density
    E_psi_e1 = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
    E_psi_e2 = 0.5*F_e.*theta_e(:,4).^2+mu^2*(1-cos(theta_e(:,3)));
    epsilon01 = E_psi_e1(end)*u_m_e(end)^1.5;
    epsilon02 = E_psi_e2(end)*u_m_e(end)^1.5;
    f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
    f2 = kappa21*f1*epsilon01/epsilon02;   %final interation for beta^2
    delta1 = 1 - E_psi_e1/4;  %delta_omega1 prepared for next stage
    delta2 = 1 - E_psi_e2/4/mu^2;  %delta_omega2 prepared for next stage
    del_psi_e1 = F_e.*theta_e(:,2).*theta_e(:,8)+sin(theta_e(:,1)).*theta_e(:,7)...
        -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,7)./u_m_e;
    del_psi_e2 = F_e.*theta_e(:,4).*theta_e(:,10)+mu^2*sin(theta_e(:,3)).*theta_e(:,9)...
        -Phi_e.*F_e.*theta_e(:,4).^2+1.5.*F_e.*theta_e(:,4).*theta_e(:,9)./u_m_e;
    Delta_theta_e1 = del_psi_e1./E_psi_e1;
    Delta_theta_e2 = del_psi_e2./E_psi_e2;
    Delta_theta_e = (f1*del_psi_e1+f2*del_psi_e2)./(f1*E_psi_e1 + f2*E_psi_e2);
    g_theta_e1 = 1.5.*F_e.*theta_e(:,2).*theta_e(:,7)./u_m_e./E_psi_e1;
    g_theta_e2 = 1.5.*F_e.*theta_e(:,4).*theta_e(:,9)./u_m_e./E_psi_e2;
    g_theta_e = (f1*E_psi_e1.*g_theta_e1+f2*E_psi_e2.*g_theta_e2)./...
        (f1*E_psi_e1+f2*E_psi_e2);
    clear del_psi_e1 del_psi_e2
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
    %%%component 2
    [pk1,loc1] = findpeaks(Delta_theta_e2.*u_m_e.^p,u_m_e);
    [pk2,loc2] = findpeaks(-Delta_theta_e2.*u_m_e.^p,u_m_e);
    index_mid1 = find(abs(loc1-u_mmid)<=2*pi/mu,1,'last');
    loc1 = loc1(1:index_mid1);
    pk1 = pk1(1:index_mid1).*loc1.^-p;
    index_mid2 = find(abs(loc2-u_mmid)<=2*pi/mu,1,'last');
    loc2 = loc2(1:index_mid2);
    pk2 = pk2(1:index_mid2).*loc2.^-p;
    pk = Delta_theta_e2(index_mid);
    diff1 = gradient(pk1)./gradient(loc1);
    diff2 = gradient(pk2)./gradient(loc2);
    %{
        We find that some low k mode would generate incorrect slope since
        findpeaks function would sometimes find both peaks and crests but 
        we want peaks only. Since gradient function subtracts nearby 2 
        points, it's better to use the penultimate point which guarantees 
        it subtracts both peaks or both crests.
    %}
    Delta_theta_e_p2 = 0.5*(diff1(end-1)-diff2(end-1));
    if A_k/mu <= 0.001
        Delta_theta_e_p2 = (pk1(end-1)-pk1(end-10001))/...
            (loc1(end-1)-loc1(end-10001));
    end
    
    Sc = 3*F_emid_sq/2/u_mmid/mu^2*(3/32/delta2(index_mid)/F_emid_sq/u_mmid^2+...
        A_k/2/F_emid_sq/u_mmid/delta2(index_mid)-...
        E_psi_e2(index_mid)/4/F_emid_sq/delta2(index_mid));
    Ss = 3/16/delta2(index_mid)/F_emid_sq/u_mmid^2/mu+...
        A_k/F_emid_sq/u_mmid/delta2(index_mid)/mu+3*F_emid_sq/2/u_mmid^2/mu-...
        3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1)/mu;
    Dp = Delta_theta_e_p2-3/2/u_mmid/delta2(index_mid)*Phi_e(index_mid)-...
        4*Phi_pe(index_mid)/delta2(index_mid);
    phis02 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid/mu*Sc);
    phic02 = 0.5*(pk+3*F_emid_sq/2/u_mmid/mu*phis02);
    
    clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p2 F_emid_sq ...
        delta2 Phi_pe
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
    g_theta_e = g_theta_e(1:index_mid);
    g_theta_e1 = g_theta_e1(1:index_mid);
    g_theta_e2 = g_theta_e2(1:index_mid);
    del_ph_p = gradient(theta_e(:,5))./gradient(u_m_e);
    del_b_p = gradient(theta_e(:,12))./gradient(u_m_e);
    %%%reshape exact solution to indexmid
    
    
    %%%ode for weak approximation
    [u_m_a,theta_a] = ode45(@(u_m,theta)...
        calculate_bytheta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
    alpha_nu,epsilon01,epsilon02,f1,f2,mu,kappa),...
        [u_mmid u_mf],...
        [theta_e(end,5);del_ph_p(end);phic01;phis01;phic02;phis02...
        ;theta_e(end,11);theta_e(end,12);del_b_p(end)]);
    clear del_ph_p del_b_p
    %%%ode for weak approximation
    
    %%%covariant energy density and velocity potential
    F_a = 1 + kappa*sqrt(u_m_a/A_m)+8/3*f1*sqrt(u_m_a)*epsilon01+...
        8/3*f2*sqrt(u_m_a)*epsilon02;
    Delta_theta_a1 = 2.*theta_a(:,3)-3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
    Delta_theta_a2 = 2.*theta_a(:,5)-3*sqrt(F_a)/2./u_m_a/mu.*theta_a(:,6);
    Delta_theta_a = (f1*epsilon01*u_m_a.^-1.5.*Delta_theta_a1+...
        f2*epsilon02*u_m_a.^-1.5.*Delta_theta_a2)./...
        (f1*epsilon01*u_m_a.^-1.5+f2*epsilon02*u_m_a.^-1.5);
    g_theta_a1 = 1.5*sqrt(F_a).*theta_a(:,4)./u_m_a;
    g_theta_a2 = 1.5*sqrt(F_a).*theta_a(:,6)./u_m_a/mu;
    g_theta_a = (f1*epsilon01*g_theta_a1+f2*epsilon02*g_theta_a2)/...
        (f1*epsilon01+f2*epsilon02);
    
    D_r = theta_a(:,2)-theta_a(:,1)./2./u_m_a;
    g_r_a = D_r./...
        (2*A_k/3./F_a+(3+1./F_a)/4./u_m_a);
    g_b_a = theta_a(:,9)./(2*A_k/3./F_a+(3+1./F_a)/4./u_m_a);
    %%%covariant energy density
    
    %%%Phi for approx
    if theta_e(index_mid,11) == 1 %if it's not turn on before mid
        index_nuturn = find(theta_a(:,7) == 1,1,'last');
        theta_off = theta_a(1:index_nuturn,:);
        u_m_off = u_m_a(1:index_nuturn,:);
        F_off = F_a(1:index_nuturn,:);
        theta_on = theta_a(index_nuturn+1:end,:);
        u_m_on = u_m_a(index_nuturn+1:end,:);
        F_on = F_a(index_nuturn+1:end,:);
        
        Phi_off = f1*u_m_off.^(-1/2).*epsilon01./A_k.*...
            (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
            f2*u_m_off.^(-1/2).*epsilon02./A_k.*...
            (2.*theta_off(:,5)-sqrt(F_off)*3/2./u_m_off/mu.*theta_off(:,6))+...
            3/8/A_k./u_m_off.*(theta_off(:,1)+kappa*theta_off(:,8).*sqrt(u_m_off/A_m));
        Phi_on = f1*u_m_on.^(-1/2)*epsilon01./A_k.*...
            (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
            f2*u_m_on.^(-1/2)*epsilon02./A_k.*...
            (2.*theta_on(:,5)-sqrt(F_on).*3/2./u_m_on/mu.*theta_on(:,6))+...
            3/8/A_k./u_m_on.*...
            (theta_on(:,1)./(1+alpha_nu)+kappa*theta_on(:,8).*sqrt(u_m_on/A_m));
        Phi_a = [Phi_off ; Phi_on];
        clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
    else
        Phi_a = f1*u_m_a.^(-1/2).*epsilon01/A_k.*...
            (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
            f2*u_m_a.^(-1/2).*epsilon02/A_k.*...
            (2.*theta_a(:,5)-sqrt(F_a)*3/2./u_m_a/mu.*theta_a(:,6))+...
            3/8/A_k./u_m_a.*...
            (theta_a(:,1)/(1+alpha_nu)+kappa*theta_a(:,8).*sqrt(u_m_a/A_m));
    end
    %%%Phi for approx
    
    %%%reshape all
    u_m = [u_m_e;u_m_a];
    Phi = [Phi_e;Phi_a];
    Delta_theta = [Delta_theta_e;Delta_theta_a];
    Delta_theta1 = [Delta_theta_e1;Delta_theta_a1];
    Delta_theta2 = [Delta_theta_e2;Delta_theta_a2];
    g_theta = [g_theta_e;g_theta_a];
    g_theta1 = [g_theta_e1;g_theta_a1];
    g_theta2 = [g_theta_e2;g_theta_a2];
    Delta_b = [theta_e(:,12);theta_a(:,8)];
    Delta_ph = [theta_e(:,5);theta_a(:,1)];
    g_b = [theta_e(:,13);g_b_a];
    g_ph = [theta_e(:,6);g_r_a];
    %%%reshape all
else
    %%%theta ode exact
    opts = odeset('MaxStep',res/mu);
    [u_m_e,theta_e] = ode45(@(u_m,theta)...
        calculate_bytheta_strong(u_m,theta,A_k,A_m,alpha_b,alpha_nu,...
        f1,f2,mu,kappa,coef_photon),[u_m0 u_mmid],...
        [Theta_i1;Theta_ip1;Theta_i2;Theta_ip2;Del_ph0;g_ph0;...
        psi01;psi11;psi02;psi12;nu_turnon],opts);
    Del_b = 3/4*theta_e(:,5);
    g_b = 3/4*theta_e(:,6);
    theta_e = [theta_e Del_b g_b];
    clear Del_b g_b
    %%%theta ode exact
    
    %%% calculate phi
    index_nuturn = find(theta_e(:,11) == 1,1,'last');
    F_e = (3/4+3/4*kappa*sqrt(u_m_e./A_m)+2*f1*u_m_e.^2.*(1-cos(theta_e(:,1)))...
        +mu^2*2*f2*u_m_e.^2.*(1-cos(theta_e(:,3))))...
        ./(0.75-f1*u_m_e.^2.*theta_e(:,2).^2-f2*u_m_e.^2.*theta_e(:,4).^2);
    theta_off = theta_e(1:index_nuturn,:);
    u_m_off = u_m_e(1:index_nuturn,:);
    F_off = F_e(1:index_nuturn,:);
    theta_on = theta_e(index_nuturn+1:end,:);
    u_m_on = u_m_e(index_nuturn+1:end,:);
    F_on = F_e(index_nuturn+1:end,:);
    
    Phi_nuoff = (u_m_off.*F_off.*theta_off(:,2).*theta_off(:,8)+...
        u_m_off.*sin(theta_off(:,1)).*theta_off(:,7)+...
        1.5.*F_off.*theta_off(:,2).*theta_off(:,7)+...
        f2/f1*(u_m_off.*F_off.*theta_off(:,4).*theta_off(:,10)+...
        mu^2*u_m_off.*sin(theta_off(:,3)).*theta_off(:,9)+...
        1.5.*F_off.*theta_off(:,4).*theta_off(:,9))+...
        +3/8/f1./u_m_off.*(theta_off(:,5)...
        +kappa*theta_off(:,12).*sqrt(u_m_off/A_m)))./...
        (u_m_off.*F_off.*(theta_off(:,2).^2+f2/f1*theta_off(:,4).^2)...
        -A_k/f1);
    Phi_nuon = (u_m_on.*F_on.*theta_on(:,2).*theta_on(:,8)...
        +u_m_on.*sin(theta_on(:,1)).*theta_on(:,7)+...
        1.5.*F_on.*theta_on(:,2).*theta_on(:,7)+...
        f2/f1*(u_m_on.*F_on.*theta_on(:,4).*theta_on(:,10)...
        +mu^2*u_m_on.*sin(theta_on(:,3)).*theta_on(:,9)+...
        1.5.*F_on.*theta_on(:,4).*theta_on(:,9))...
        +3/8/f1./u_m_on.*(theta_on(:,5)/(1+alpha_nu)...
        +kappa*theta_on(:,12).*sqrt(u_m_on/A_m)))...
        ./(u_m_on.*F_on.*(theta_on(:,2).^2+f2/f1*theta_on(:,4).^2)...
        -A_k/f1);
    Phi_e = [Phi_nuoff ; Phi_nuon];
    Phi_pe = gradient(Phi_e)./gradient(u_m_e);
    clear theta_off u_m_off F_off theta_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
    %%% calculate phi
    
    %%%covariant energy density
    E_psi_e1 = 0.5*F_e.*theta_e(:,2).^2+1-cos(theta_e(:,1));
    E_psi_e2 = 0.5*F_e.*theta_e(:,4).^2+mu^2*(1-cos(theta_e(:,3)));
    epsilon01 = E_psi_e1(end)*u_m_e(end)^1.5;
    epsilon02 = E_psi_e2(end)*u_m_e(end)^1.5;
    f1 = 3*(1-kappa)/8/sqrt(u_meq)/epsilon01/(1+kappa21);
    f2 = kappa21*f1*epsilon01/epsilon02;   %final interation for beta^2
    delta1 = 1 - E_psi_e1/4;  %delta_omega1 prepared for next stage
    delta2 = 1 - E_psi_e2/4/mu^2;  %delta_omega2 prepared for next stage
    del_psi_e1 = F_e.*theta_e(:,2).*theta_e(:,8)+sin(theta_e(:,1)).*theta_e(:,7)...
        -Phi_e.*F_e.*theta_e(:,2).^2+1.5.*F_e.*theta_e(:,2).*theta_e(:,7)./u_m_e;
    del_psi_e2 = F_e.*theta_e(:,4).*theta_e(:,10)+mu^2*sin(theta_e(:,3)).*theta_e(:,9)...
        -Phi_e.*F_e.*theta_e(:,4).^2+1.5.*F_e.*theta_e(:,4).*theta_e(:,9)./u_m_e;
    Delta_theta_e1 = del_psi_e1./E_psi_e1;
    Delta_theta_e2 = del_psi_e2./E_psi_e2;
    Delta_theta_e = (f1*del_psi_e1+f2*del_psi_e2)./(f1*E_psi_e1 + f2*E_psi_e2);
    g_theta_e1 = 1.5.*F_e.*theta_e(:,2).*theta_e(:,7)./u_m_e./E_psi_e1;
    g_theta_e2 = 1.5.*F_e.*theta_e(:,4).*theta_e(:,9)./u_m_e./E_psi_e2;
    g_theta_e = (f1*E_psi_e1.*g_theta_e1+f2*E_psi_e2.*g_theta_e2)./...
        (f1*E_psi_e1+f2*E_psi_e2);
    clear del_psi_e1 del_psi_e2
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
    
    clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p1 delta1 
    clear Sc Ss Dp
    %%%component 2
    [pk1,loc1] = findpeaks(Delta_theta_e2.*u_m_e.^p,u_m_e);
    [pk2,loc2] = findpeaks(-Delta_theta_e2.*u_m_e.^p,u_m_e);
    index_mid1 = find(abs(loc1-u_mmid)<=2*pi/mu,1,'last');
    loc1 = loc1(1:index_mid1);
    pk1 = pk1(1:index_mid1).*loc1.^-p;
    index_mid2 = find(abs(loc2-u_mmid)<=2*pi/mu,1,'last');
    loc2 = loc2(1:index_mid2);
    pk2 = pk2(1:index_mid2).*loc2.^-p;
    pk = Delta_theta_e2(index_mid);
    diff1 = gradient(pk1)./gradient(loc1);
    diff2 = gradient(pk2)./gradient(loc2);
    %{
        We find that some low k mode would generate incorrect slope since
        findpeaks function would sometimes find both peaks and crests but 
        we want peaks only. Since gradient function subtracts nearby 2 
        points, it's better to use the penultimate point which guarantees 
        it subtracts both peaks or both crests.
    %}
    Delta_theta_e_p2 = 0.5*(diff1(end-1)-diff2(end-1));
    if A_k/mu <= 0.001
        Delta_theta_e_p2 = (pk1(end-1)-pk1(end-10001))/...
            (loc1(end-1)-loc1(end-10001));
    end
    
    
    Sc = 3*F_emid_sq/2/u_mmid/mu^2*(3/32/delta2(index_mid)/F_emid_sq/u_mmid^2+...
        A_k/2/F_emid_sq/u_mmid/delta2(index_mid)-...
        E_psi_e2(index_mid)/4/F_emid_sq/delta2(index_mid));
    Ss = 3/16/delta2(index_mid)/F_emid_sq/u_mmid^2/mu+...
        A_k/F_emid_sq/u_mmid/delta2(index_mid)/mu+3*F_emid_sq/2/u_mmid^2/mu-...
        3/8/u_mmid^2/F_emid_sq*(F_e(index_mid)-1)/mu;
    Dp = Delta_theta_e_p2-3/2/u_mmid/delta2(index_mid)*Phi_e(index_mid)-...
        4*Phi_pe(index_mid)/delta2(index_mid);
    phis02 = (2*Dp-Sc*pk)/(2*Ss+3*F_emid_sq/2/u_mmid/mu*Sc);
    phic02 = 0.5*(pk+3*F_emid_sq/2/u_mmid/mu*phis02);
    
    %clear pk pk1 pk2 loc1 loc2 diff1 diff2 Delta_theta_e_p2 F_emid_sq ...
    %    delta2 Phi_pe index_mid1 index_mid2
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
    g_theta_e = g_theta_e(1:index_mid);
    g_theta_e1 = g_theta_e1(1:index_mid);
    g_theta_e2 = g_theta_e2(1:index_mid);
    del_ph_p = gradient(theta_e(:,5))./gradient(u_m_e);
    del_b_p = gradient(theta_e(:,12))./gradient(u_m_e);
    %%%reshape exact solution to indexmid
    
    
    %%%ode for weak approximation
    [u_m_a,theta_a] = ode45(@(u_m,theta)...
        calculate_bytheta_weak_appox(u_m,theta,A_k,A_m,alpha0,beta0,...
    alpha_nu,epsilon01,epsilon02,f1,f2,mu,kappa),...
        [u_mmid u_mtr],...
        [theta_e(end,5);del_ph_p(end);phic01;phis01;phic02;phis02...
        ;theta_e(end,11);theta_e(end,12);del_b_p(end)]);
    clear del_ph_p del_b_p
    %%%ode for weak approximation
    
    %%%covariant energy density
    F_a = 1 + kappa*sqrt(u_m_a/A_m)+8/3*f1*sqrt(u_m_a)*epsilon01+...
        8/3*f2*sqrt(u_m_a)*epsilon02;
    Delta_theta_a1 = 2.*theta_a(:,3)-3*sqrt(F_a)/2./u_m_a.*theta_a(:,4);
    Delta_theta_a2 = 2.*theta_a(:,5)-3*sqrt(F_a)/2./u_m_a/mu.*theta_a(:,6);
    Delta_theta_a = (f1*epsilon01*u_m_a.^-1.5.*Delta_theta_a1+...
        f2*epsilon02*u_m_a.^-1.5.*Delta_theta_a2)./...
        (f1*epsilon01*u_m_a.^-1.5+f2*epsilon02*u_m_a.^-1.5);
    %%%covariant energy density
    
    %%%Phi for approx
    if theta_e(index_mid,11) == 1 %if it's not turn on before mid
        index_nuturn = find(theta_a(:,7) == 1,1,'last');
        theta_off = theta_a(1:index_nuturn,:);
        u_m_off = u_m_a(1:index_nuturn,:);
        F_off = F_a(1:index_nuturn,:);
        theta_on = theta_a(index_nuturn+1:end,:);
        u_m_on = u_m_a(index_nuturn+1:end,:);
        F_on = F_a(index_nuturn+1:end,:);
        
        Phi_off = f1*u_m_off.^(-1/2).*epsilon01./A_k.*...
            (2.*theta_off(:,3)-sqrt(F_off)*3/2./u_m_off.*theta_off(:,4))+...
            f2*u_m_off.^(-1/2).*epsilon02./A_k.*...
            (2.*theta_off(:,5)-sqrt(F_off)*3/2./u_m_off/mu.*theta_off(:,6))+...
            3/8/A_k./u_m_off.*(theta_off(:,1)+kappa*theta_off(:,8).*sqrt(u_m_off/A_m));
        Phi_on = f1*u_m_on.^(-1/2)*epsilon01./A_k.*...
            (2.*theta_on(:,3)-sqrt(F_on).*3/2./u_m_on.*theta_on(:,4))+...
            f2*u_m_on.^(-1/2)*epsilon02./A_k.*...
            (2.*theta_on(:,5)-sqrt(F_on).*3/2./u_m_on/mu.*theta_on(:,6))+...
            3/8/A_k./u_m_on.*...
            (theta_on(:,1)./(1+alpha_nu)+kappa*theta_on(:,8).*sqrt(u_m_on/A_m));
        Phi_a = [Phi_off ; Phi_on];
        clear psi_off u_m_off F_off psi_on u_m_on F_on Phi_nuoff Phi_nuon index_nuturn
    else
        Phi_a = f1*u_m_a.^(-1/2).*epsilon01/A_k.*...
            (2.*theta_a(:,3)-sqrt(F_a)*3/2./u_m_a.*theta_a(:,4))+...
            f2*u_m_a.^(-1/2).*epsilon02/A_k.*...
            (2.*theta_a(:,5)-sqrt(F_a)*3/2./u_m_a/mu.*theta_a(:,6))+...
            3/8/A_k./u_m_a.*...
            (theta_a(:,1)/(1+alpha_nu)+kappa*theta_a(:,8).*sqrt(u_m_a/A_m));
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
    
    %%%2nd phase velocity
    g_theta_a1 = 1.5*sqrt(F_a).*theta_a(:,4)./u_m_a;
    g_theta_a2 = 1.5*sqrt(F_a).*theta_a(:,6)./u_m_a/mu;
    g_theta_a = (f1*epsilon01*g_theta_a1+f2*epsilon02*g_theta_a2)/...
        (f1*epsilon01+f2*epsilon02);
    
    D_r = theta_a(:,2)-theta_a(:,1)./2./u_m_a;
    g_r_a = D_r./...
        (2*A_k/3./F_a+(3+1./F_a)/4./u_m_a);
    g_b_a = theta_a(:,9)./(2*A_k/3./F_a+(3+1./F_a)/4./u_m_a);
    %%%2nd phase velocity
    
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
    g_r0 = (D_r0-A_tr*u_mtr^-1.5*sqrt(F_a(end))*Q0)/...
        (2*A_k/3/F_a(end)+(3+1/F_a(end))/4/u_mtr);
    %%%initial condition for 3rd phase
    
    %%%ode 3rd phase
    [u_m_3,theta_3] = ode45(...
    @(u_m,theta) calculate_bytheta_weak_re(u_m,theta,A_k,A_m,A_t,alpha_b,alpha_nu...
    ,epsilon01,epsilon02,f1,f2,mu,m_e,BE,T_CMB,u_mn,n_b0,kappa),...
        [u_mtr u_mf],...
        [theta_a(end,1);g_r0;theta_a(end,3);theta_a(end,4)...
        ;theta_a(end,5);theta_a(end,6);theta_a(end,7);theta_a(end,8);Q0]);
    %%%ode 3rd phase
    
    %%%some parameters
    F_3 = 1 + kappa*sqrt(u_m_3/A_m)+8./3.*sqrt(u_m_3)*f1*epsilon01...
        +8./3.*sqrt(u_m_3)*f2*epsilon02;
    Phi_3 = -epsilon01*f1/A_k.*u_m_3.^(-1/2).*(2.*theta_3(:,3)...
        -3.*sqrt(F_3)/2./u_m_3.*theta_3(:,4))...
        -epsilon02*f2/A_k.*u_m_3.^(-1/2).*(2.*theta_3(:,5)...
        -3.*sqrt(F_3)/2./u_m_3/mu.*theta_3(:,6))...
        -3/8/A_k./u_m_3.*(theta_3(:,1)/(1+alpha_nu)+...
        kappa*sqrt(u_m_3/A_m).*theta_3(:,8));
    %%%some parameters
    
    %%%covariant energy density/velocity 3
    Delta_theta_31 = 2.*theta_3(:,3)-sqrt(F_3)*3/2./u_m_3.*theta_3(:,4);
    Delta_theta_32 = 2.*theta_3(:,5)-sqrt(F_3)*3/2./u_m_3/mu.*theta_3(:,6);
    Delta_theta_3 = (f1*epsilon01*u_m_3.^-1.5.*Delta_theta_31+...
        f2*epsilon02*u_m_3.^-1.5.*Delta_theta_32)./...
        (f1*epsilon01*u_m_3.^-1.5+f2*epsilon02*u_m_3.^-1.5);
    g_theta_31 = 1.5*sqrt(F_3).*theta_3(:,4)./u_m_3;
    g_theta_32 = 1.5*sqrt(F_3).*theta_3(:,6)./u_m_3/mu;
    g_theta_3 = (f1*epsilon01*g_theta_31+f2*epsilon02*g_theta_32)/...
        (f1*epsilon01+f2*epsilon02);
    %%%covariant energy density/velocity 3
     
    %%%reshape all
    u_m = [u_m_e;u_m_a;u_m_3];
    Phi = [Phi_e;Phi_a;Phi_3];
    Delta_theta = [Delta_theta_e;Delta_theta_a;Delta_theta_3];
    Delta_theta1 = [Delta_theta_e1;Delta_theta_a1;Delta_theta_31];
    Delta_theta2 = [Delta_theta_e2;Delta_theta_a2;Delta_theta_32];
    Delta_b = [theta_e(:,12);theta_a(:,8);theta_3(:,8)];
    Delta_ph = [theta_e(:,5);theta_a(:,1);theta_3(:,1)];
    g_b = [theta_e(:,13);g_b_a;0.75*(theta_3(:,2)+theta_3(:,9))];
    g_ph = [theta_e(:,6);g_r_a;theta_3(:,2)];
    g_theta = [g_theta_e;g_theta_a;g_theta_3];
    g_theta1 = [g_theta_e1;g_theta_a1;g_theta_31];
    g_theta2 = [g_theta_e2;g_theta_a2;g_theta_32];
    %%%reshape all
    
end
toc

%%%plot
figure(1)
loglog(sqrt(u_m/u_meq),abs(Delta_theta))
hold on
loglog(sqrt(u_m/u_meq),abs(Delta_b))
hold on
loglog(sqrt(u_m/u_meq),abs(Delta_ph))
legend('\PsiDM_{total}','Baryon_\Psi','Photon_\Psi','location','Southeast');
title('coariant energy density')
xlabel('a/a_{eq}')
ylabel('\Delta')
set(gca,'FontSize',16)

figure(2)
loglog(sqrt(u_m/u_meq),abs(g_theta))
hold on
loglog(sqrt(u_m/u_meq),abs(g_b))
hold on
loglog(sqrt(u_m/u_meq),abs(g_ph))
legend('\PsiDM_{total}','Baryon_\Psi','Photon_\Psi','location','Southeast');
title('velocity potential')
xlabel('a/a_{eq}')
ylabel('g')
set(gca,'FontSize',16)

figure(3)
loglog(sqrt(u_m/u_meq),abs(Delta_theta))
hold on
loglog(sqrt(u_m/u_meq),abs(Delta_theta1))
hold on
loglog(sqrt(u_m/u_meq),abs(Delta_theta2))
legend('\PsiDM_{total}','\PsiDM_1','\PsiDM_2','location','Southeast');
xlabel('a/a_{eq}')
ylabel('\Delta')
title('component density')
set(gca,'FontSize',16)
xlabel('a/a_{eq}')

figure(4)
loglog(sqrt(u_m/u_meq),abs(g_theta))
hold on
loglog(sqrt(u_m/u_meq),abs(g_theta1))
hold on
loglog(sqrt(u_m/u_meq),abs(g_theta2))
legend('\PsiDM_{total}','\PsiDM_1','\PsiDM_2','location','Southeast');
xlabel('a/a_{eq}')
ylabel('g')
title('component velocity')
set(gca,'FontSize',16)


figure(5)
loglog(sqrt(u_m/u_meq),abs(Phi))
hold on
xlabel('a/a_{eq}')
ylabel('\phi')
title('gravity')
set(gca,'FontSize',16)


%{
file = fopen('EA_appox.txt','wt');
fprintf(file,'%24.17f \t %24.17f \n',[u_m,abs(Delta_theta)]');
fclose(file);
%}
%{
legend('\PsiDM','Baryon_\Psi','Photon_\Psi','CDM','Baryon_{CDM}','Photon_{CDM}','location','Southeast')
xlabel('a/a_{eq}')
title('covariant energy density \Theta_0 = 179^\circ , k =0.1k_c')

legend('Axion','CDM')
xlabel('a/a_{eq}')
title('Potential \Theta_0 = 179^\circ , k =0.1k_c')
%}
%%%plot
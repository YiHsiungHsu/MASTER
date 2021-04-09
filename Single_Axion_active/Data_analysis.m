clear all
close all

%%%input parameter
%End time in a_eq (scale factor at matter-radiation equality)
a_f = 30;

%Initial wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_initial = 0.1;

%End wave mode in critical k_c22 (typically k_c22~(100 kpc)^-1)
k_end = 2;

%particle mass in 10^-22 eV
mr = 1;

%background axion angle shift from the top of the potential hill:
angle_shift = 90;

%number of headerlines
HeaderLines_CDM = 7;
HeaderLines_FDM = 8;

%%%input parameter

%%%load file
afstr = num2str(a_f);
ki = num2str(k_initial);
kf = num2str(k_end);
m = num2str(mr*1e-22);
angle = num2str(180 - angle_shift);

%%%CDM
filename = ['CDM_' afstr 'aeq' '_from_k_' ki 'to' kf '.txt'];
fid = fopen(filename,'r');
formatSpec = '%f %f %f %f %f %f %f';
CDM = textscan(fid,formatSpec,'HeaderLines',HeaderLines_CDM);
CDM = cell2mat(CDM);
%%%CDM

%%%FDM
filename = ['Axion_' m '_' angle '_' afstr 'aeq' '_from_k_' ki 'to' kf '.txt'];
fid = fopen(filename,'r');
formatSpec = '%f %f %f %f %f %f %f';
FDM = textscan(fid,formatSpec,'HeaderLines',HeaderLines_FDM);
FDM = cell2mat(FDM);
%%%FDM

k = FDM(:,1);
TF = (FDM./CDM).^2;

DM_transfer = TF(:,2);
baryon_transfer = TF(:,3);
photon_transfer =TF(:,4);
g_DM_transfer = TF(:,5);
g_baryon_transfer = TF(:,6);
g_photon_transfer =TF(:,7);

figure(1)
semilogy(k,DM_transfer,'Color',[0, 0, 1],'linewidth',2)
xmin = k_initial;
xmax = k_end;
xlim([xmin,xmax])
grid on
set(gca,'FontSize',14)

figure(2)
semilogy(k,g_DM_transfer,'Color',[0, 0, 1],'linewidth',2)
xmin = k_initial;
xmax = k_end;
xlim([xmin,xmax])
grid on
set(gca,'FontSize',14)
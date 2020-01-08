clear
close all
clc
delete(gcp('nocreate'))

%% Declare constants
C = Constants;
%% Load functions
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/othercolor';
addpath(fDIR)

%% data
skala=100;
teta=0;
R=(0.5:0.5:20)*10^(-6);
delta=0.00001:0.0005:0.01001;

%%

[RR,deltas]=meshgrid(R,delta);
tau_ps=(2*C.ro_p*RR.^2)./(9*C.nu*C.ro_a);
z_b=C.g*tau_ps*cos(teta).*deltas.^2/C.nu;

figHandler = figure('Color',[1 1 1]);
width=16;
height=20;
fsize=13;
set(gcf,'Position', [300, 300, 2*560, 2.3*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

surf(deltas*skala,RR*10^6,z_b*skala)
xlabel('\delta[cm]')
ylabel('R [um]');
zlabel('z_b [cm]')
xlim([0 1.01001])
%ylim([])
set(gca,'FontSize',fsize)
%set(gca,'clim',[1, 40])
box on
grid on
h=colorbar
ylabel(h, 'z_b [cm]')
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
R=(0.5:0.5:20)*10^(-6);
%delta=0.00001:0.0005:0.01001;
delta=0.001:0.0005:0.021;

%%

[RR,deltas]=meshgrid(R,delta);
tau_ps=(2*C.ro_p*RR.^2)./(9*C.nu*C.ro_a);
k=(sqrt(1+(4*C.nu*tau_ps./(deltas.^2)))).^(-1);
tau_z=2*tau_ps./(k.^(-1)-1);

figHandler = figure('Color',[1 1 1]);
width=16;
height=20;
fsize=13;
set(gcf,'Position', [300, 300, 2*560, 2.3*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

surf(deltas*skala,RR*10^6,tau_z)
xlabel('\delta [cm]')
ylabel('$R [\mu m]$','interpreter','latex');
zlabel('\tau_z [s]')
%xlim([])
%ylim([])
set(gca,'FontSize',fsize)
%set(gca,'clim',[1, 40])
box on
grid on
h=colorbar
ylabel(h, '\tau_z [s]')
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
teta=pi/4;
R=(1:0.01:25)*10^(-6);
delta=0.001:0.00005:0.01001;

%%

[RR,deltas]=meshgrid(R,delta);
tau_ps=(2*C.ro_p*RR.^2)./(9*C.nu*C.ro_a);
z_b=C.g*tau_ps*cos(teta).*deltas.^2/C.nu;

figHandler = figure('Color',[1 1 1]);
width=16;
height=20;
fsize=16;
set(gcf,'Position', [300, 300, 2*560, 2.3*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off',...
    'Color','white')
colormap jet
pcolor(deltas*skala,RR*10^6,z_b*skala)
xlabel('$\delta$~[cm]','interpreter','latex','FontSize',fsize+2)
ylabel('$R$~[$\mu$m]','interpreter','latex','FontSize',fsize+2);
zlabel('$z_b$~[cm]','interpreter','latex','FontSize',fsize+2)
xlim([0.1 1.0001])
%ylim([])
set(gca,'FontSize',fsize)
%set(gca,'clim',[1, 40])
box on
shading interp
h=colorbar
ylabel(h,'$z_b$~[cm]','interpreter','latex','FontSize',fsize+2)
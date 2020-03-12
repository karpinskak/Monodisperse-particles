clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';

%% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
% plot parameters
fsize=18;
width=16;
height=20;

%% parameter ranges
R=(1:20)*10^(-6);
delta=(0.001:0.01:1)*2.5*10^(-2);
[RR,Delta]=meshgrid(R,delta);
teta=0.5*pi/2;

tau_p=2*Const.ro_p.*RR.^2/(9*Const.nu*Const.ro_a);

P=1+Const.nu^(-1)*Const.g^2*tau_p.^3*sin(teta)^2-Const.nu*tau_p.*Delta.^(-2);
C=(4*pi)^(-2);
Amax=(0.5*C*P.*((1+4*Const.nu*tau_p./(P.^2.*Delta.^2)).^(1/2)-1)).^(1/2);
Amax=real(Amax).*(imag(Amax)==0);
Amax(Amax==0)=NaN;

figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'Color',[1 1 1],...
    'papersize',[width,height],...
    'InvertHardCopy','off')

pcolor(Delta*10^2,RR*10^6,Amax)
pp=colorbar;
ylabel(pp,'$A_{max}$','interpreter','latex','FontSize',fsize)
shading interp
hold on
contour(Delta*10^2,RR*10^6,Amax,[Const.Acr,Const.Acr],'--','ShowText','on','Color','k')
text(2.152,19,'$\theta=\frac{\pi}{4}$','interpreter','latex','FontSize',fsize)

ylabel('$R [\mu m]$','interpreter','latex')
xlabel('$\delta [cm]$','interpreter','latex')
ylim([1 20])
set(gca,'FontSize',fsize)
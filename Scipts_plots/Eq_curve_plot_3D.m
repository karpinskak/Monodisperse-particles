clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';

%% Load constants and functions
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
ri=2.1866;
% give A values for plot
AA=unique([10^(-8),0.0000001:0.00000001:0.000002]);
rplot=[0.00001,0.001:4:600];
% plot parameters
fsize=16;
width=16;
height=20;

%% plot
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')
[rr,Aa]=meshgrid(rplot,AA);
%eq_curve(rr,Aa)

surface(rr,Aa,eq_curve(rr,Aa))
%shading interp

%ylim([0 0.2])
xlabel('$r^{\ast}$','interpreter','latex')
ylabel('$A$','interpreter','latex')
zlabel('$f_A(r^{\ast})$','interpreter','latex')
%ticzki=get(gca, 'XTick');
%ticz=unique([Const.rs, ri, ticzki]);
%set(gca, 'XTick', ticz,...
%'FontSize',fsize);
%legend()
hold off
grid on
box on

function krzywa=eq_curve(r,A)

chi=(1-exp(-r.^2/2))./(2*pi*A.*r.^2);
krzywa=A.*r.*sqrt(1+chi.^2);
end
%%
clear
close all
clc
delete(gcp('nocreate'))
tic

DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];

% plot parameters
fsize=16;
width=16;
height=20;
skala=100;

%% Declare constants, load data
Const = Constants;
name='Dock_time_eps_00001_full.mat';
load([loadDIR name])
addpath([DIR 'Functions/'])
R=13*10^(-6);
[par]=wylicz_param(Const,3,R,ST,AA,0,0,0);
tau_f=par.tau_f;

%% Grid data
X=AA;
Y=ST;
Z=reshape(texitt,size(AA)).*tau_f;%.*par.gamma/2;

%% plot
close all
figure
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

pcolor(X,Y,Z)
 colormap(jet)
 c=colorbar;
 set(gca,'ColorScale','log')
set(gca,'clim',[min(min(Z)) max(max(Z))])
shading interp
xlabel('A')
ylabel('St')
ylim([0 1])
xlim([min(A) max(A)])
c.Label.String='t_{doc} [s]';
set(gca,'FontSize',fsize)
hold off
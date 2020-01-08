clear
close all
clc
delete(gcp('nocreate'))
tic

%% Declare constants, load functions and data
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';

loadDIR=[DIR spec];
Const = Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
load([loadDIR 'Dock_time_eps_00001.mat'])
%%
figure(1)
tdocest=-log(eps/(Const.rs-eps))./Aa;
y=Stt;
z1=real(tdocest);
z2=real(texit2);
yplot=y;
z1plot=z1;
z2plot=z2;

% surf(Aa.*wybiorcze,yplot,z1plot)
% xlabel('A')
% ylabel('St/A')
% zlabel('t_{doc}-est')
% figure(2)
% surf(Aa.*wybiorcze,yplot,z2plot)
% xlabel('A')
% ylabel('St/A')
% zlabel('t_{doc}')
% zlim(8*10^4)

surf(Aa,yplot,log(abs(z2plot-z1plot)./z1plot))
xlabel('A')
 ylabel('St/A')
 zlabel('t_{num}-t_{est}')
 colorbar
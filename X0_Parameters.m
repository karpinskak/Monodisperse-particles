% choose directory
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/';
addpath(fDIR)
name='Dock_time_eps_00001_inorbit_detTraj_testkatow.mat';

%% Zasoby
poolnr=12;

% Simulation parameters
% give A and St as a list of paired parameters
%AA=0.0001:0.001:0.0101;
%ST=[0.0001:0.0001:0.001,0.0011:0.001:0.101,0.111:0.01:1];
%AA=[0.0001:0.0001:0.0081];
ST=0.0001:0.05:1.0001;
katy=linspace(atan(160),0.9995*pi/2,25);
katy(end)=[];
St_A=tan(katy);
AA=ST'./St_A;
Stt=repmat(ST',1,numel(St_A));

A=reshape(AA,numel(AA),1);
St=reshape(Stt,numel(AA),1);


t00=0;
eps=0.0001;
nt=2000;
lim_time=6000; %[s]
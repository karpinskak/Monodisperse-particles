% choose directory
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/';
addpath(fDIR)
name='Dock_time_eps_00001_detTraj.mat';

%% Zasoby
poolnr=13;

% Simulation parameters
% give A and St as a list of paired parameters
%AA=0.0001:0.001:0.0101;
%ST=[0.0001:0.0001:0.001,0.0011:0.001:0.101,0.111:0.01:1];
%AA=[0.0001:0.0001:0.0081];
ST=0.0001:0.05:1.0001;
St_A=160:500:9900;
AA=ST'./St_A;
Stt=repmat(ST',1,numel(St_A));

A=reshape(AA,numel(AA),1);
St=reshape(Stt,numel(AA),1);


t00=0;
eps=0.0001;
nt=1000;
tfinn=5000; % initial 'final simulation time' - changed later if not sufficient
lim_time=3600; %[s]
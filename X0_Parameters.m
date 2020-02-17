% choose directory
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/';
addpath(fDIR)
name='Dock_time_eps_00001_full.mat';

%% Zasoby
poolnr=13;

% Simulation parameters
% give A and St as a list of paired parameters
%AA=0.0001:0.001:0.0101;
st=[0.0001:0.02:1];
aa=[0.0001:0.0001:0.0065];
[ST,AA]=meshgrid(st,aa);
A=reshape(AA,numel(AA),1);
St=reshape(ST,numel(AA),1);

% full i katy
% st1=0.0001:0.05:1.0001;
% katy=[atan(10):(atan(158)-atan(10))/25:atan(158),atan(160):(atan(160)-atan(10))/25:atan(1000)];
% St_A=tan(katy);
% AA1=st1'./St_A;
% STT1=repmat(st1',1,numel(St_A));
% A1=reshape(AA1,numel(AA1),1);
% St1=reshape(STT1,numel(AA1),1);
% St1(A1>0.01)=[];
% A1(A1>0.01)=[];
% st2=0.0011:0.01:1;
% aa2=0.0001:0.0005:0.01;
% [STT2,AA2]=meshgrid(st2,aa2);
% A2=reshape(AA2,numel(AA2),1);
% St2=reshape(STT2,numel(AA2),1);
% A=[A1;A2];
% St=[St1;St2];


t00=0;
eps=0.0001;
nt=2000;
lim_time=7200; %[s]
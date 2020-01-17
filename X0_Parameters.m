% choose directory
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/';
addpath(fDIR)
name='Dock_time_eps_00001_full_i_katy.mat';

%% Zasoby
poolnr=14;

% Simulation parameters
% give A and St as a list of paired parameters
%AA=0.0001:0.001:0.0101;
%ST=[0.0001:0.0001:0.001,0.0011:0.001:0.101,0.111:0.01:1];
%AA=[0.0001:0.0001:0.0081];
st1=0.0001:0.05:1.0001;
katy=[atan(10):(atan(158)-atan(10))/25:atan(158),atan(160):(atan(160)-atan(10))/25:atan(1000)];
St_A=tan(katy);
AA1=st1'./St_A;
STT1=repmat(st1',1,numel(St_A));
A1=reshape(AA1,numel(AA1),1);
St1=reshape(STT1,numel(AA1),1);
St1(A1>0.01)=[];
A1(A1>0.01)=[];

st2=0.0011:0.01:1;
aa2=0.0001:0.0005:0.01;
[AA2,STT2]=meshgrid(st2,aa2);
A2=reshape(AA2,numel(AA2),1);
St2=reshape(STT2,numel(AA2),1);

A=[A1;A2];
St=[St1;St2];


t00=0;
eps=0.0001;
nt=2000;
lim_time=6000; %[s]
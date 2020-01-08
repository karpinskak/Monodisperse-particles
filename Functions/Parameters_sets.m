function [DIR,spec,poolnr,tstart,T,par_set,a,b,c,d,dispersity,Rmin,Rdev,PD,type,k,l,n,max_R,new_sim,loadPosDIR,loadDIR] = Parameters(p,sets)
%PARAMETERS Provides parameters for Simulation_control
%  This is the same as in Parameters.m but in the form of the function that
%  takes the parameters from the table.

% choose directory
DIR =sets(p).DIR;
spec=sets(p).spec;
%spec='Results/A00076/Delta0012/R12_9/Th05/';%025 dla Th oznacza 0.25*pi/2
%% Zasoby
poolnr=sets(p).poolnr;

%% Time
tstart=sets(p).tstart; %[s]
T=sets(p).T; % czas trwania symulacji [s]


%% Vortex and particle parameters
% choose input parameter set:
% 0 - for tracers
% 1 - standard nondimnesional numbers {St,Sv,theta,A}
% 2 - dimensional {R, delta, theta,A} - this case can be polydispersed

par_set=sets(p).par_set;

% Simulation set 1
a=sets(p).a; %St or R [m];
b=sets(p).b; % Sv or delta [m]
c=sets(p).c; % theta [rad]
d=sets(p).d; %A [1]

dispersity=sets(p).dispersity;
% 0- monodispersed
% 1 - uniform up to Rmean+2*Rdev,
% 2 - polydisperse gaussian
% 3- polydisperse other

switch dispersity
    case 0
        Rdev=[];
        Rmin=[];
        PD=[];
    case 1
        Rdev=sets(p).Rdev;
        Rmin=sets(p).Rmin;
        PD=[];
    case 2
        Rdev=sets(p).Rdev;
        Rmin=sets(p).Rmin;
        PD(:,1) =str2num(sets(p).PD);
        PD(:,2) =normpdf(PD(:,1),a,Rdev)*(PD(2,1)-PD(1,1));
    case 3
        load('/home/pracownicy/karpinska/Dokumenty/Ruch_kropli/Dziury_dane/Dane_od_Tiny/Probab_dist_27.mat','PD')
        PD(:,1)=PD(:,1)*10^(-6);
        a=sum(PD(:,1).*PD(:,2));
        Rmin=sets(p).Rmin;
        Rdev=[];
end

%Small test
if dispersity~=0 && par_set~=2
    error('There is no combination of this parameter set with polidispersity.')
end

%% Domain filling
type=sets(p).type; % 0 - randomly generated positions, initially full cylindrical domain
% 1 - evenly generatel positions at time t=0,

k=sets(p).k; % D=k*delta, radius
l=sets(p).l; % Z=l*delta, [-Z,Z]
n=sets(p).n; %[nr of particles/m^3]
max_R=sets(p).max_R; % particle initial positions random generator grid size scale [m]

%% saved and loaded data
new_sim=sets(p).new_sim; % 0 - generate new initial particle positions, 1- use existing
% initial particle positions, 2 - restart old simulation
% (inactive)

if type==0
    nam='';
elseif type==1
    nam='_even';
end

if par_set==0
    nam2='_tracer';
else
    nam2='';
end

switch dispersity
    case 2
        nam3='_poli_gauss';
    case 3
        nam3='_poli_exper';
    case 1
        nam3='_poli_uni';
    case 0
        nam3='_mono';
end
% Format k number to put it in the folder name
kt = round(k);
k1=sprintf('%d',kt);
if abs(kt-k)<10^(-6)
   str = sprintf(['k', k1]);
else
    k2t = (k-round(k));
    k2=sprintf('%0.2f',k2t);
    k2=k2(3:end);
    str = sprintf(['k', k1, '_', k2]);
end

loadPosDIR=[DIR spec str '_l' num2str(l) '_n' num2str(n/1000000) '_tk' num2str(tstart+T) nam nam3];
loadDIR=[loadPosDIR nam2];
end


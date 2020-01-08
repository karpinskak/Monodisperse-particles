
% choose directory
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results/A00012/St025/Sv0016/Th05/';
%025 dla Th oznacza 0.25*pi/2
%% Zasoby
poolnr=10;

%% Time
tstart=0; %[s]
T=4; % czas trwania symulacji [s]


%% Vortex and particle parameters
% choose input parameter set:
% 0 - for tracers
% 1 - standard nondimnesional numbers {St,Sv,theta,A}
% 2 - dimensional {R, delta, theta,A} - this case can be polydispersed

par_set=1;

% Simulation set 1
 a=0.25; %St or R [m];
 b=0.016; % Sv or delta [m]
 c=0.5; % theta [rad]
 d=0.0012; %A [1]


dispersity=0;
% 0- monodispersed
% 1 - uniform up to Rmax,
% 2 - polydisperse gaussian
% 3- polydisperse other

switch dispersity
    case 1
        Rmax=12*10^(-6);
        Rmin=2*10^(-6);
    case 2
        Rdev=4.85*10^(-6);
        Rmin=1*10^(-6);
        PD(:,1) =(1:0.5:90)*10^(-6);
        PD(:,2) = normpdf(PD(:,1),a,Rdev)*(PD(2,1)-PD(1,1));
    case 0
        Rdev=0;
        Rmin=0;
    case 3
        load('/home/pracownicy/karpinska/Dokumenty/Ruch_kropli/Dziury_dane/Dane_od_Tiny/Probab_dist_27.mat','PD')
        PD(:,1)=PD(:,1)*10^(-6);
        a=sum(PD(:,1).*PD(:,2));
        Rmin=1*10^(-6);
end

%Small test
if dispersity~=0 && par_set~=2
    error('There is no combination of this parameter set with polidispersity.')
end

%% Domain
type=1; % 0 - randomly generated positions, initially full cylindrical domain
% 1 - evenly generatel positions at time t=0,

k=6; % D=k*delta, radius
l=2; % Z=l*delta, [-Z,Z]
n=0.1*10^6; %[nr of particles/m^3]
max_R=100*10^(-6); % particle initial positions random generator grid size scale [m]

%% saved and loaded data
new_sim=0; % 0 - generate new initial particle positions, 1- use existing
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
clear nam nam2 nam3

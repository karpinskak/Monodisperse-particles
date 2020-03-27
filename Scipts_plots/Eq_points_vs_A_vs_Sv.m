clc
delete(gcp('nocreate'))
tic

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';

%% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)

% variable ranges
AA=0.00001:2*10^(-5):0.03001;
Svv=unique([0.01:0.01:0.1,0.072,0.0815]);

% plot parameters
% xmax=0.023;
% ymax=28;
% al=0.5;
fsize=16;
width=16;
height=20;
% style={'-','--','-'};
% colors=hsv(numel(Svv));

%% calc
equ_points=cell([numel(Svv),1]);
parpool('local',3)
parfor j=1:numel(Svv)
    Sv_temp=Svv(j);
    punkty=zeros(numel(AA),6);
    for k=1:numel(AA)
        A_temp=AA(k);
        point_temp = eq_points(A_temp,Sv_temp);
        [np,~]=size(point_temp);
        punkty(k,1:2*np)=reshape(point_temp',[1,2*np]);
        if np<3
            punkty(k,(2*np+1):6)=NaN;
        end
        
    end
    equ_points{j}=punkty;
end

%% reshape eq. points
punkt1=zeros(size(Svv,A));
punkt2=zeros(size(Svv,A));
punkt3=zeros(size(Svv,A));
[~,n_Acr]=min(abs(AA-Const.Acr));
if AA(n_Acr)>Const.Acr
    n_Acr=n_Acr-1;
end
[Sv,A]=meshgrid(Svv,AA);

for k=1:n_Acr
    for j=1:numel(Svv)
        punkty=equ_points{j};
    end
end
for k=(n_Acr+1):numel(AA)
    for j=1:numel(Svv)
        
    end
end
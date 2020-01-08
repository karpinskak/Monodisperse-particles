clear
close all
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
AA=unique([0.001:0.002:0.031,0.02176]);
%AA=unique([0.001:0.005:0.031,0.02176]);
Svv=10^(-4):10^-(3):0.9001;
%Svv=10^(-4):10^-(3):0.9001;
St=0.5;

% plot parameters
xmax=0.023;
ymax=28;
al=0.5;
fsize=16;
width=16;
height=20;
style={'-','--','-'};
colors=jet(numel(Svv));

%% calc
equ_points=cell([numel(AA),1]);
parpool('local',3)
parfor j=1:numel(AA)
    A_temp=AA(j);
    punkty=zeros(numel(Svv),6);
    for k=1:numel(Svv)
        Sv_temp=Svv(k);
        point_temp = eq_points(A_temp,Sv_temp);
        [np,~]=size(point_temp);
        punkty(k,1:2*np)=reshape(point_temp',[1,2*np]);
        if np<3
            punkty(k,(2*np+1):6)=NaN;
        end
        
    end
    equ_points{j}=punkty;
end

toc

%% plot
clf
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

leg=cell([numel(Svv),1]);

for j=1:numel(AA)
    A_temp=AA(j);
    leg{j}=num2str(A_temp);
    dane=equ_points{j};
    prom=(dane(:,1:2:6).^2+dane(:,2:2:6).^2).^(1/2);
    for k=1:3
        if k==1
            plot(Svv,prom(:,k),style{k},'Color',colors(j,:),'LineWidth',2)
        else
            plot(Svv,prom(:,k),style{k},'Color',colors(j,:),'LineWidth',2,'HandleVisibility','off')
        end
        hold on
    end
%     % plot unstable points as crosses in the domain r<r_s
%     I=find((prom<Const.rs));
%     AAi=AA(I);
%     promi=prom(I);
%     Stcr=AAi./abs(fun_fi(promi))';
%     J=find(Stcr<St);
%     if j==numel(Svv)
%         scatter(AAi(J(1:10:end)),promi(J(1:10:end)),10,'x','MarkerFaceColor','none','MarkerEdgeColor',[0 0 0]);
%     leg{j+1}='unstable';
%     else
%     scatter(AAi(J(1:10:end)),promi(J(1:10:end)),10,'x','MarkerFaceColor','none','MarkerEdgeColor',[0 0 0],'HandleVisibility','off');
%     end
%     hold on
%     clear I J promi AAi
end
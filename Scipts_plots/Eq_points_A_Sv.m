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
AA=0.00001:2*10^(-5):0.03001;
Svv=unique([0.01:0.01:0.09,0.072,0.0815]);
St=0.5;

% plot parameters
xmax=0.023;
ymax=28;
al=0.5;
fsize=16;
width=16;
height=20;
style={'-','--','-'};
colors=hsv(numel(Svv));

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

toc

%% plot
clf
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

% plot stability regions
rplot=Const.rs:0.001:ymax;
fiplot=fun_fi(rplot).^(1/2);
patch([0,xmax,xmax,0],[0,0,Const.rs,Const.rs],[0.72 0.27 0.1],'EdgeColor','none','FaceAlpha',al-0.1,'HandleVisibility','off')
patch([fiplot,xmax,xmax,0],[rplot,ymax,Const.rs,Const.rs],[0.6 0.65 0.7],'EdgeColor','none','FaceAlpha',al-0.2,'HandleVisibility','off')
patch([fiplot,0,0],[rplot,ymax,Const.rs],[0.6 0.7 0.7],'EdgeColor','none','FaceAlpha',al,'HandleVisibility','off')
hold on
annotation('textbox',[0.56 0.66 0.2 0.2],'String','stable node','Color',[0.55 0.6 0.65],'LineStyle','none','FontSize',fsize+2,'FitBoxToText','on');
annotation('textbox',[0.13,0.22,0.10,0.046],'String','saddle','Color',[0.55 0.65 0.65],'LineStyle','none','FontSize',fsize+2,'FitBoxToText','on');
annotation('textbox',[0.67,0.12,0.26,0.046],'String','stable/unstable focus','Color',[0.75 0.3 0.13],'LineStyle','none','FontSize',fsize+2,'FitBoxToText','on');


leg=cell([numel(Svv),1]);

% plot equilibrium points r_0(A)
for j=1:numel(Svv)
    Sv_temp=Svv(j);
    leg{j}=num2str(Sv_temp);
    dane=equ_points{j};
    prom=(dane(:,1:2:6).^2+dane(:,2:2:6).^2).^(1/2);
    for k=1:3
        if k==1
            plot(AA,prom(:,k),style{k},'Color',colors(j,:),'LineWidth',2)
        else
            plot(AA,prom(:,k),style{k},'Color',colors(j,:),'LineWidth',2,'HandleVisibility','off')
        end
        hold on
    end
    % plot unstable points as crosses in the domain r<r_s
    I=find((prom<Const.rs));
    AAi=AA(I);
    promi=prom(I);
    Stcr=AAi./abs(fun_fi(promi))';
    J=find(Stcr<St);
    if j==numel(Svv)
        scatter(AAi(J(1:10:end)),promi(J(1:10:end)),10,'x','MarkerFaceColor','none','MarkerEdgeColor',[0 0 0]);
    leg{j+1}='unstable';
    else
    scatter(AAi(J(1:10:end)),promi(J(1:10:end)),10,'x','MarkerFaceColor','none','MarkerEdgeColor',[0 0 0],'HandleVisibility','off');
    end
    hold on
    clear I J promi AAi
end

% formatting
ylim([0 ymax])
xlim([0 xmax])
  
ylabel('$r^{\ast}_0$','interpreter','latex')
xlabel('A','interpreter','latex')
set(gca,'FontSize',fsize+2)

a=legend(leg);
 a.Title.String='S_v';
 a.Title.FontSize=fsize-2;
 a.FontSize=fsize-4;
 hold on

hold off
grid on

ticzkiX=get(gca, 'XTick');
ticzX=unique([Const.Acr, ticzkiX]);
I=find(ticzX==Const.Acr);
set(gca, 'XTick', ticzX,...
    'FontSize',fsize);
tlabX=get(gca, 'XTickLabel');
tlabX{I}='A_{cr}';
set(gca,'XTickLabel', tlabX)

ticzkiY=get(gca, 'YTick');
ticzY=unique([Const.rs, Const.ri, ticzkiY]);
J=find(ticzY==Const.rs);
K=find(ticzY==Const.ri);
set(gca, 'YTick', ticzY,...
    'FontSize',fsize);
tlabY=get(gca, 'YTickLabel');
tlabY{J}='r_s';
tlabY{K}='r_i';
set(gca,'YTickLabel', tlabY)

clear tlabX tlab tlabY ticzX ticzY ticzkiX ticzkiY I J K

%% functions
function wart=fun_fi(r)
part1=(1-exp(-r.^2/2))./r.^2;
part3=exp(-r.^2/2);
wart=part1.*(part1-part3)/(2*pi)^2;
end
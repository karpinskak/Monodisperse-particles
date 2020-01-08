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

npool=14;

%% define scenarios [Sv,A]
Sv=[.06,.06,.06,.06,.06,.06,.078,.078,.0815];
A=[.008,.01132,.013,.008,.01132,.013,.0189,.0142,.0228];
St=[0.001,0.001,0.001,5,5,5,.5,.5,.5];
 
%% give plot parameters
D=5;
ticz=-6:2:6;
N=5;
tspan=0:0.1:500;

% figure parameters
fsize=16;
width=16;
height=20;
p_style={'o','o','x','s'};
opt={'none','k','none','k'};
opt2={'k','none','k','none'};
tytuly={'1a','2a','3a','1b','2b','3b','4','5','6'};

%% calc trajectory
[pol_r,pol_fi]=positions_even(D,N);
M1=2*A./St;
M2=St.^(-1);
M3=Sv./St;
parpool('local',npool)
traj=cell(numel(A),1);
parfor k=1:numel(A)
    traj_temp=cell(numel(pol_r),1);
    for p=1:numel(pol_r)
        trajektorie=zeros(numel(tspan),3);
        x0=[pol_r(p),pol_fi(p),0,0];
        [t1,x1]=ode45(@(t,x)trajektoria_polarne(t,x,M1(k),M2(k),M3(k)),tspan,x0);
    trajektorie(:,1)=x1(:,1).*cos(x1(:,2));
    trajektorie(:,2)=x1(:,1).*sin(x1(:,2));
    trajektorie(:,3)=tspan;
    traj_temp{p}=trajektorie;
    end
    traj{k}=traj_temp;
    disp(['Zestaw ' num2str(k) 'sie policzyl.'])
end
delete(gcp('nocreate'))

%% stability points
punkty=zeros(numel(A),6);
stability=zeros(numel(A),3);
for j=1:numel(A)
    [point_temp, stab_temp] = eq_points(Const,A(j),Sv(j),St(j));
    [np,~]=size(point_temp);
    punkty(j,1:2*np)=reshape(point_temp',[1,2*np]);
    stability(j,1:np)=stab_temp;
    if np<3
        punkty(j,(2*np+1):6)=NaN;
        stability(j,(np+1):3)=0;
    end
end

toc
save('scenarios_9p_T_500_D_5_N_5.mat','-v7.3')
toc
%% plot
clf
kolory=jet(numel(pol_r));
figure(1)
set(gcf,'Position', [534,1,1022,970],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')


tiledlayout(3, 3)

for  k=1:numel(A)
    traj_temp=traj{k};
    ax(k)=nexttile;
    for n=1:numel(pol_r)
        temp=traj_temp{n};
    plot(temp(:,1),temp(:,2),'Color',kolory(n,:))
    hold on
    end
    for n=1:3
        scatter(punkty(k,2*n-1),punkty(k,2*n),p_style{stability(k,n)+1},'LineWidth',2,'MarkerFaceColor',opt{stability(k,n)+1},'MarkerEdgeColor',opt2{stability(k,n)+1})
    end
    
    %title(['A=' num2str(A(k)) ', Sv=' num2str(Sv(k))])
    title(['(' tytuly{k} ')'])
    xlabel('x^+')
    ylabel('y^+')
    axis equal
    grid on
    set(gca, 'XTick', ticz);
    hold off
end
axis(ax,[-6 6 -8 4])
tlt.Padding = "none";
tlt.TileSpacing = "none";

 %% functions
 function [pol_r,pol_fi]=positions_even(D,N)
ile_kaw=N-1;
Xlin=-D:2*D/ile_kaw:D;
Ylin=-D:2*D/ile_kaw:D;
[X,Y]=meshgrid(Xlin,Ylin);
R=(X.^2+Y.^2).^(1/2);
FI=atan2(Y,X);
pol_r=reshape(R,[numel(Xlin)*numel(Ylin),1]);
pol_fi=reshape(FI,[numel(Xlin)*numel(Ylin),1]);
if numel(pol_r)~=numel(pol_fi)
    error('Cos nie tak z losowaniem jednorodnym')
end
 end

% Approximate radial velocity with exponent for axis docking process

clear
close all
clc
delete(gcp('nocreate'))
tic

%% Declare constants, load functions and data
DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';

loadDIR=[DIR spec];
Const = Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)

load([loadDIR 'Dock_time_eps_00001_2.mat'])

%% Calc and plot parameters
poolnr=12;
options=optimset('MaxFunEvals',100000000);
% plot parameters
fsize=16;
width=16;
height=20;

%% Approximation calculation

n=3;
typ=1;
switch typ
    case 1
        r0=1.58499;
    case 0
        r0=0;
end

wybor=[];
for p=1:numel(part)
    param(p)=part(p).par.A^(-1)*part(p).par.St;
    if abs(part(p).init.r0_bez-r0)<=eps && imag(texit(p))==0 && isnan(texit(p))==0
        wybor=[wybor,p];
    end
end

param_wyb=param(wybor);
parmax=max(param_wyb);

parpool('local',poolnr)

parfor k=1:numel(wybor)
    
    xplot=part(wybor(k)).traj.t;
    yplot= part(wybor(k)).traj.X(:,n);
    % cutting off some false points?
    
    [rad_vel_field,~,~]=velocity_field(part(wybor(k)).par.A,Const.rs-eps,0);
    
    yplot= yplot./rad_vel_field;

    
    % yplot=a*exp(-xplot/b)+c
    fun = @(x)sseval3(x,xplot,yplot);
    x0=[1,part(wybor(k)).par.St/part(wybor(k)).par.A,0];
    [wyn,~,exitflag] = fminsearch(fun,x0,options);
    if exitflag==0
        fun_param(k)=NaN;
        mnoz_param(k)=NaN;
        wolny_param(k)=NaN;
        Rsq(k)=NaN;
    else
        fun_param(k)=wyn(2);
        mnoz_param(k)=wyn(1);
        wolny_param(k)=wyn(3);
        yplot2=wyn(1)*exp(-xplot/wyn(2))+wyn(3);
        Rsq(k) = 1 - sum((yplot - yplot2).^2)/sum((yplot2 - mean(yplot2)).^2);
    end
    
end
toc
delete(gcp('nocreate'))
save('tau_d1_fit_eps_00001_2.mat')
    
%% plot

set(0,'DefaultAxesColor','none')
clf
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
            'paperunits','centimeters',...
            'papersize',[width,height],...
            'InvertHardCopy','off')
        
do_plot=(Rsq<=1).*(Rsq>0);
Awyb=A(wybor)';
Stwyb=St(wybor)';
Yplot=fun_param(do_plot==1).*Awyb(do_plot==1);
Xplot=(param_wyb(do_plot==1));%.^(-1);
sizes=floor(5^2*Awyb(do_plot==1)/max(AA));
scatter(Xplot,Yplot,sizes,[zeros(size(Stwyb(do_plot==1))); Stwyb(do_plot==1)/max(Stwyb(do_plot==1)); 1- Stwyb(do_plot==1)/max(Stwyb(do_plot==1)) ]','o','filled')
set(gca,'FontSize',fsize)
xlabel('St/A','interpreter','latex')
ylabel('$\tau^+_{d 1}*A$','interpreter','latex')
box on
grid on
xlim([0 160])
ylim([1,1.25])
cb=colorbar
title(cb,'St')
%%
function sse = sseval3(x,tdata,ydata)
a = x(1);
b = x(2);
c= x(3);
sse = sum((ydata - (a*exp(-tdata/b)+c)).^2);
end
function sse=sseval_new(x,tdata,ydata)
b=x(1);
c=x(2);
d=x(3);
m=x(4);
n=x(5);
sse = sum((ydata - (b*(1-exp(-(x-c)^n/d))/(tdata-c)^m)).^2);
end
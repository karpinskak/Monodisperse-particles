clear
close all
clc
delete(gcp('nocreate'))
tic

DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/Initial_vel_fluid/';
loadDIR=[DIR spec];

%% Declare constants, load functions and data
Const = Constants;
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/';
addpath(fDIR)

load([loadDIR 'Dock_time_eps_00001.mat'])

%% choose plot type
typ=0;  % 0- z eps do r_orb
        % 1- z r_s do eps,

% how many St/A cases - plots in one figure
ile=15;
%
fsize=16;
width=20;
height=25;

%skalowanie: 0- nie, 1- tak
skal=0;

%
switch typ
    case 1
        r0=1.58499;
    case 0
        r0=0;
end

%% plot
% choose relevant St/A values (for respective docking process)
wybor=[];
param=zeros(1,numel(part));
for p=1:numel(part)
    param(p)=part(p).par.A^(-1)*part(p).par.St;
    if abs(part(p).init.r0_bez-r0)<=eps && imag(texit(p))==0 && isnan(texit(p))==0
        wybor=[wybor,p];
    end
end
% choose just "ile" values for plotting
param_wyb=param(wybor);
parmin=min(param_wyb);
parmax=max(param_wyb);
wartosci=[];
wartosci=linspace(parmin,parmax,ile);

if typ==1 && skal==0
wartosci(1)=[];
wartosci(end)=[];
end
numerki=[];
for j=1:numel(wartosci)
    [out,idx]=sort(abs(param_wyb-wartosci(j)));
   znalezione=find(param_wyb(idx(1))==param);
    numerki(j)=znalezione(1);
end
numerki=unique(numerki,'stable');
if typ==0
numerki(1)=[];
end

param_rys=param(numerki);
maxt=max(texit(numerki));

% plot
figure
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')
 set(gca,'FontSize',fsize)
 
for k=1:numel(numerki)
    p=numerki(k);
    leg{k}=num2str(round(param_rys(k)));
    switch typ
        case 0
            switch skal
                case 1
                    xplot=part(p).traj.t/texit(p);
                    labelkax='$t^{\ast}/t_{doc}$';
                    limitx=[0 1];
                case 0
                    xplot=part(p).traj.t;
                    labelkax='$t^{\ast}$';
                    limitx=[0 maxt];
            end
            
            subplot(2,2,1)
            
            %rorb=Row_orb_Bur(param(p));
            %yplot1=part(p).traj.X(:,1)/rorb;
            yplot1=part(p).traj.X(:,1);
            plot(xplot,yplot1,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            ylabel('$r^{\ast}$','interpreter','latex')
            text(-0.2,1.1,'a','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
            set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,2)
            yplot2=part(p).traj.X(:,2)/(2*pi);
            plot(xplot,yplot2,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            ylabel('$\phi^{\ast}/2 \pi$','interpreter','latex')
            text(-0.2,1.1,'b','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,3)
            yplot3= part(p).traj.X(:,3);
            [rad_vel_field,~,~]=velocity_field(part(p).par.A,eps,0);
            if skal==1
                yplot3= yplot3./rad_vel_field;
            end
            plot(xplot,yplot3,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            if skal==0
                ylabel('$\dot{r^{\ast}}$','interpreter','latex')
            else
                ylabel('$\dot{r^{\ast}}/v_r^{\ast}(\epsilon)$','interpreter','latex')
            end
            text(-0.2,1.1,'c','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,4)
            yplot4= part(p).traj.X(:,4);
            if skal==1
                yplot4=yplot4/sqrt(param(p)^(-1));
            end
            plot(xplot,yplot4,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            if skal==0
                ylabel('$\dot{\phi^{\ast}}$','interpreter','latex')
            else
                ylabel('$\dot{\phi^{\ast}}/\omega_{orb}$','interpreter','latex')
            end
            text(-0.2,1.1,'d','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
        case 1
            
            %  n==3
            % [rad_vel_field,~,~]=velocity_field(part(p).par.A,Const.rs,0);
            % yplot= -xplot./log(yplot./rad_vel_field);
            
            switch skal
                case 1
                    xplot=part(p).traj.t/texit(p);
                    labelkax='$t^{\ast}/t_{doc}$';
                    limitx=[0 1];
                case 0
                    xplot=part(p).traj.t;
                    labelkax='$t^{\ast}$';
                    limitx=[0 maxt];
            end
            
            subplot(2,2,1)
            yplot1=part(p).traj.X(:,1);
            plot(xplot,yplot1,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            ylabel('$r^{\ast}$','interpreter','latex')
            text(-0.2,1.1,'a','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,2)
            yplot2=part(p).traj.X(:,2)/(2*pi);
            plot(xplot,yplot2,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            ylabel('$\phi^{\ast}/2 \pi$','interpreter','latex')
            text(-0.2,1.1,'b','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,3)
            yplot3= part(p).traj.X(:,3);
            [rad_vel_field,~,~]=velocity_field(part(p).par.A,Const.rs,0);
            if skal==1
            yplot3= yplot3./rad_vel_field;
            end
            plot(xplot,yplot3,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            if skal==0
                ylabel('$\dot{r^{\ast}}$','interpreter','latex')
            else
                ylabel('$\dot{r^{\ast}}/v_r^{\ast}$','interpreter','latex')
            end
            text(-0.2,1.1,'c','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
                 set(gca,'FontSize',fsize)
            hold on
            
            subplot(2,2,4)
            yplot4= part(p).traj.X(:,4);
            if skal==1
             [~,az_vel_field,~]=velocity_field(part(p).par.A,eps,0);
            yplot4=yplot4/az_vel_field;
            end
            plot(xplot,yplot4,'Color',[param(p)/parmax 1-param(p)/parmax 0])
            xlabel(labelkax,'interpreter','latex');
            xlim(limitx)
            if skal==0
                ylabel('$\dot{\phi^{\ast}}$','interpreter','latex')
            else
            ylabel('$\dot{\phi^{\ast}}/\frac{v_{\phi}(\epsilon)}{\epsilon}$','interpreter','latex')
            end
            text(-0.2,1.1,'d','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            grid on
            set(gca,'FontSize',fsize)
            hold on
    end
    legenda=legend(leg)
    title(legenda,'St/A')
    set(gca,'FontSize',fsize)%,'Position',[1.2,-2.2, 0.0786, 0.456])
end
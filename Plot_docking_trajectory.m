%%% A script that plots a 4-tiled figure with r,dor{r}, fi, dot{fi}, scaled
%%% or not. Tracking data for multiple parameters from group "numerki" is
%%% ploted in the same figure. Used by W_Part_tracking_plots_corr.mlx

% choose labels and axis limits for all the tiles
switch skal
    case 1
        labelkax='$t^+/t_{doc}^+$';
        limitx=[0 1];
        switch typ
            case 0
                labelkay='$\dot{r^+}/v_r^+(\epsilon)$';
                labelkay2='$\dot{\phi^+}/\omega_{orb}$';
            case 1
                labelkay='$\dot{r^+}/v_r^+$';
                labelkay2='$\dot{\phi^+}/\frac{v_{\phi}(\epsilon)}{\epsilon}$';
        end
    case 0
        labelkax='$t^+$';
        labelkay='$\dot{r^+}$';
        labelkay2='$\dot{\phi^+}$';
        limitx=[0 maxt];
end
% colorscale
kolory=[kol_odn/max(kol_odn) 1-kol_odn/max(kol_odn) zeros(size(kol_odn))];


%%%%%% plot %%%%%%
tiledlayout(4, 1)
% 1. tile
nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,yplot1,~,~,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot1,'Color',kolory(p,:))
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel('$r^+$','interpreter','latex')

% 2. tile
nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,yplot2,~,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot2,'Color',kolory(p,:))
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel('$\phi/2 \pi$','interpreter','latex')

% 3. tile
nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,~,yplot3,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot3,'Color',kolory(p,:))
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay,'interpreter','latex')

% 4. tile
nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,~,~,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot4,'Color',kolory(p,:))
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay2,'interpreter','latex')
   

set(gca,'FontSize',fsize)

% choose relevant data for the plot and scale it if needed
function [xplot,yplot1,yplot2,yplot3,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p)
 xplot=part(p).traj.t;
    yplot1=part(p).traj.X(:,1);
    yplot2=part(p).traj.X(:,2)/(2*pi);
    yplot3= part(p).traj.X(:,3);
    yplot4= part(p).traj.X(:,4);
       
    if skal==1
        xplot=xplot/texit(p);
        switch typ
            case 0 % in-orbit
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,eps,0);
                yplot3= yplot3./rad_vel_field;
                [~,az_vel_field,~]=velocity_field(part(p).par.A,yplot1,0);
                yplot4=yplot4./az_vel_field;
                %yplot4=yplot4/sqrt(param(p)^(-1));

            case 1 % point
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,Const.rs,0);
                yplot3= yplot3./rad_vel_field;
                [~,az_vel_field,~]=velocity_field(part(p).par.A,eps,0);
                yplot4=yplot4/az_vel_field;
                
        end
    end
end

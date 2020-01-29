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
                labelkay0='$r^+/r^+_{orb}$';
                labelkay1='$-\dot{r^+}/u_r^+(\sigma)$';
                labelkay2='$r^+\varphi(r^+)/u^+_{\varphi}(r)$';
                %labelkay2='$\dot{\phi^+}/\omega_{orb}$';
            case 1
                labelkay0='$r^+$';
                labelkay1='$\dot{r^+}/u_r^+$';
                labelkay2='$\sigma\dot{\varphi^+}/u^+_{\varphi}(\sigma)$';
        end
    case 0
        labelkax='$t^+$';
        labelkay0='$r^+$';
        labelkay1='$\dot{r^+}$';
        labelkay2='$\dot{\varphi^+}$';
        limitx=[0 maxt];
end
%%%% colorscales %%%%%

%convert basic color from rgb to hsl
colb_hsl=rgb2hsl(col_basic(p,:));
kolory_hsl=[colb_hsl(1)+zeros(numel(numerki),1) colb_hsl(2)+zeros(numel(numerki),1) linspace(0.95, 0.3, numel(numerki))'];
%convert back
kolory=hsl2rgb(kolory_hsl);



%%%%%% plot %%%%%%

% 1. tile
nexttile(rysunek,1)
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,yplot1,~,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plocik=plot(xplot,yplot1,'Color',kolory(l,:),'Linewidth',2);
    if l==round(numel(numerki)/2)
        subsecik=[subsecik,plocik];
    end
    hold on
end
xlim(limitx)
xticklabels([])
ylabel(labelkay0,'interpreter','latex')
set(gca,'FontSize',fsize)
grid on

% 2. tile
nexttile(rysunek,2)
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,yplot3,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot3,'Color',kolory(l,:),'Linewidth',2)
    hold on
end
xlim(limitx)
xticklabels([])
%ylim(limity)
ylabel(labelkay1,'interpreter','latex')
set(gca,'FontSize',fsize)
grid on

% 3. tile
nexttile(rysunek,3)
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,~,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot4,'Color',kolory(l,:),'Linewidth',2)
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay2,'interpreter','latex')
if skal==1 && typ==0
    ylim([0.99 1.001])
end
grid on   

set(gca,'FontSize',fsize)

% choose relevant data for the plot and scale it if needed
function [xplot,yplot1,yplot3,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p)
 xplot=part(p).traj.t;
    yplot1=part(p).traj.X(:,1);
    yplot3= part(p).traj.X(:,3);
    yplot4= part(p).traj.X(:,4);
       
    if skal==1
        xplot=xplot/texit(p);
        switch typ
            case 0 % in-orbit
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,eps,0);
                yplot3= -yplot3./rad_vel_field;
                [~,az_vel_field,~]=velocity_field(part(p).par.A,yplot1,0);
                yplot4=smoothdata(yplot4./az_vel_field,'gaussian',500);
                %yplot4=yplot4/sqrt(param(p)^(-1));
                yplot1=yplot1/(max(yplot1));

            case 1 % point
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,Const.rs,0);
                yplot3= yplot3./rad_vel_field;
                [~,az_vel_field,~]=velocity_field(part(p).par.A,eps,0);
                yplot4=yplot4/az_vel_field;
                
        end
    end
end

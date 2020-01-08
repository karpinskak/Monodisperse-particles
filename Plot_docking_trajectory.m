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



tiledlayout(4, 1)

nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,yplot1,~,~,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot1,'Color',[param(p)/parmax 1-param(p)/parmax 0])
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel('$r^+$','interpreter','latex')

nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,yplot2,~,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot2,'Color',[param(p)/parmax 1-param(p)/parmax 0])
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel('$\phi/2 \pi$','interpreter','latex')


nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,~,yplot3,~]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot3,'Color',[param(p)/parmax 1-param(p)/parmax 0])
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay,'interpreter','latex')


nexttile
for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,~,~,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p);
    plot(xplot,yplot4,'Color',[param(p)/parmax 1-param(p)/parmax 0])
    hold on
end
xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay2,'interpreter','latex')
   

set(gca,'FontSize',fsize)

function [xplot,yplot1,yplot2,yplot3,yplot4]=plot_variables(skal,typ,eps,part,param,texit,p)
 xplot=part(p).traj.t;
    yplot1=part(p).traj.X(:,1);
    yplot2=part(p).traj.X(:,2)/(2*pi);
    yplot3= part(p).traj.X(:,3);
    yplot4= part(p).traj.X(:,4);
       
    if skal==1
        xplot=xplot/texit(p);
        switch typ
            case 0
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,eps,0);
                yplot3= yplot3./rad_vel_field;
                yplot4=yplot4/sqrt(param(p)^(-1));
            case 1
                [rad_vel_field,~,~]=velocity_field(part(p).par.A,Const.rs,0);
                yplot3= yplot3./rad_vel_field;
                [~,az_vel_field,~]=velocity_field(part(p).par.A,eps,0);
                yplot4=yplot4/az_vel_field;
                
        end
    end
end

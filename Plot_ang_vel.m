labelkax='$t^+/t_{doc}^+$';
limitx=[0 1];
labelkay2='$r^+\varphi(r^+)/u^+_{\varphi}(r)$';

%%%% colorscales %%%%%

%convert basic color from rgb to hsl
colb_hsl=rgb2hsl(col_basic(p,:));
kolory_hsl=[colb_hsl(1)+zeros(numel(numerki),1) colb_hsl(2)+zeros(numel(numerki),1) linspace(0.95, 0.3, numel(numerki))'];
%convert back
kolory=hsl2rgb(kolory_hsl);

set(gcf,'Position', [640, 300, 2*560, 2*420/3 ])


for l=1:numel(numerki)
    p=numerki(l);
    [xplot,~,yplot4]=plot_variables(skal,typ,part,texit,p);
    plocik=plot(xplot,yplot4,'Color',kolory(l,:),'Linewidth',2);
    hold on
    if l==round(numel(numerki)/2)
        subsecik=[subsecik,plocik];
    end
end

xlabel(labelkax,'interpreter','latex');
xlim(limitx)
ylabel(labelkay2,'interpreter','latex')
grid on   
set(gca,'FontSize',fsize)

function [xplot,yplot1,yplot4]=plot_variables(skal,typ,part,texit,p)

xplot=part(p).traj.t;
yplot1=part(p).traj.X(:,1);
yplot4= part(p).traj.X(:,4);

if skal==1
    xplot=xplot/texit(p);
    switch typ
        case 0 % in-orbit
            [~,az_vel_field,~]=velocity_field(part(p).par.A,yplot1,0);
            yplot4=smoothdata(yplot4./az_vel_field,'gaussian',500);
            
        case 1 % point
            [~,az_vel_field,~]=velocity_field(part(p).par.A, yplot1, 0);
            yplot4=smoothdata(yplot4./az_vel_field,'gaussian',250);
            
    end
end
end
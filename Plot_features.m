set(gcf,'Position', [640, 300, 1.4*560, 1.7*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')
axis equal
set(gca,'XLim',[xmin  xmax], 'YLim',[ymin ymax],...
    'FontSize',fsize,'Color','white')
xlabel('x[cm]')
ylabel('y[cm]')
if dimension==2
    set(gca,'XTick',linspace(-part(1).par.D*skala,part(1).par.D*skala,11));
set(gca,'YTick',linspace(-part(1).par.D*skala,part(1).par.D*skala,11));

elseif dimension==3
set(gca,'ZLim', [zmin zmax])
zlabel('z[cm]')
set(gca,'XTick',linspace(ceil(xmin),floor(xmax),5));
set(gca,'YTick',linspace(ceil(ymin),floor(ymax),5));
set(gca,'ZTick',linspace(ceil(zmin),floor(zmax),12));

end
box on
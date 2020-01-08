%%
clear
close all
clc
delete(gcp('nocreate'))
tic

DIR ='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';
spec='Results_dock_time/';
loadDIR=[DIR spec];

%% Declare constants
Const = Constants;

% Plot parameters
fsize=16;
width=12;
height=30;

%% load and reshape data
close all
figure
set(gcf,'Position', [640, 300,  2*400,5*280 ],...
        'paperunits','centimeters',...
        'papersize',[width,height],...
        'InvertHardCopy','off')
for k=1:5
    
    name=['Dock_time_eps_',sprintf('%0*.f',k+1,1),'.mat'];
    load([loadDIR name])
    A2=reshape(A,size(Aa));
    St2=reshape(St,size(Aa));
    texit2=reshape(texit,size(Aa));
    texit2(texit2==0)=NaN;
    odl=reshape(abs(r0_bez-rfin_bez),size(Aa));
    
    %% plot
    
    subplot(5,2,2*k-1)
    hold on
    surf(A2(:,1:30),St2(:,1:30),log(real(texit2(:,1:30))),'EdgeColor','none')
    shading interp
    view([-30 40])
    hold off
    zlim([5 15])
    caxis manual
    caxis([5 15])

    subplot(5,2,2*k)
    pcolor(A2(:,1:30),St2(:,1:30),log(real(texit2(:,1:30))))
    hold on
    AA(16*pi^2*AA>1)=[];
    plot(AA,16*pi^2*AA,'Color',[0 0 0],'LineWidth',1)
    shading interp
    caxis manual
    caxis([5 15])
    colorbar
    hold off
end
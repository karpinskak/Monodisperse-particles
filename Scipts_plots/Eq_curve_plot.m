clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';

%% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR 'Functions/'];
addpath(fDIR)
ri=2.1866;
% give A values for plot
AA=unique([0.001:0.002:0.04,0.02176]);
rplot=[0.00001,0.001:0.0001:25];
% plot parameters
fsize=16;
width=16;
height=20;

%% plot
figure(1)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

% equillibrium curves
for k=1:numel(AA)
    krzywa=@(r)eq_curve(r,AA(k));
    krzywa_plot=krzywa(rplot);
    if AA(k)==0.02176
        gl=2;
    else
        gl=1;
    end
    plot(rplot,krzywa_plot,'Color',[AA(k)/max(AA) 0 1-AA(k)/max(AA)],'Linewidth',gl)
    leg{k}=num2str(AA(k));
    hold on
end
ylim([0 0.12])
xlabel('$r^{\ast}$','interpreter','latex')
ylabel('$f_A(r^{\ast})$','interpreter','latex')

ticzki=get(gca, 'XTick');
ticz=unique([Const.rs, ri, ticzki]);
set(gca, 'XTick', ticz,...
    'FontSize',fsize);
legenda=legend(leg);
title(legenda,'A')
grid on
box on
% krzywa ograniczajaca wszystkie krzywe z dolu
plot(rplot,(1-exp(-rplot.^2/2))./(2*pi*rplot),'--','Color',[0.5 0.5 0.5])
hold off

figure(2)
set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

subplot(2,1,1)
plot(rplot,stab_curve(rplot),'LineWidth',2)
ticzki=get(gca, 'XTick');
ticz=unique([Const.rs, ri, ticzki]);
set(gca, 'XTick', ticz,...
    'FontSize',fsize);
xlabel('$r^{\ast}$','interpreter','latex')
ylabel('$\varphi(r)$','interpreter','latex')
grid on
box on
text(-0.2,1.1,'a','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(2,1,2)
yplot=sqrt(stab_curve(rplot));
yplot(real(yplot)==0)=NaN;
plot(rplot,yplot,'LineWidth',2)
ticzki=get(gca, 'XTick');
ticz=unique([Const.rs, ri, ticzki]);
set(gca, 'XTick', ticz,...
    'FontSize',fsize);
xlabel('$r^{\ast}$','interpreter','latex')
ylabel('$\sqrt{\varphi(r)}$','interpreter','latex')
grid on
box on
text(-0.2,1.1,'b','FontSize',fsize,'Units', 'Normalized', 'VerticalAlignment', 'Top')

function krzywa=eq_curve(r,A)

chi=(1-exp(-r.^2/2))./(2*pi*A*r.^2);
krzywa=A*r.*sqrt(1+chi.^2);
end

function phi=stab_curve(r)

part1=(1-exp(-r.^2/2))./r.^2;
part3=exp(-r.^2/2);
phi=part1.*(part1-part3)/(2*pi)^2;
end
% This script plots 3 elements of vortex velocity field

clear
close all
clc
delete(gcp('nocreate'))
fsize=16;
width=16;
height=20;

r_bez=0.00001:0.01:15.00001;
z_bez=0.00001:0.01:15.00001;
A=0.001;

vr_bez=-A*r_bez;
vfi_bez=(1-exp(-r_bez.^2/2))./(2*pi*r_bez.^2);
vz_bez=2*A*z_bez;

set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')


yyaxis left
plot(r_bez,-vr_bez,r_bez,vfi_bez,'LineWidth',3)
ylabel('$u^+$','interpreter','latex')
yyaxis right
plot(z_bez,vz_bez,'LineWidth',3)
ylabel('$u^+_z$','interpreter','latex')
xlabel('$r^+$','interpreter','latex')
set(gca,'FontSize',fsize)
text(12,0.032,['A=', num2str(A)],'FontSize',fsize+4)
text(4,0.018,'$-u_{\varphi}^+$','FontSize',fsize+4,'interpreter','latex','Color',[0,0.45, 0.74])
text(4.3,0.0028,'$u^+_r$','FontSize',fsize+4,'interpreter','latex','Color',[0,0.45, 0.74])
grid on
xlim([0,15])
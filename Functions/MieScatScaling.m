function [Iscat,Ix]=MieScatScaling(prom,promx,kat)
%% Mie scattering scaling
% prom is the array of particle radii [um]
% promx is the size to which we scale the rest according to Mie scattering
% intensity [um]
% kat is arbitrary Mie scattering angle
% Iscat and Ix are vectors of relative scattering intensity values for prom
% and promx respectively

addpath('/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Mie_scattering/')
load('scattering_d_theta_01_R_1_01_90.mat')
[~,I]=min(abs(scatt.angles-kat));
prom=[prom,promx];
klaster_nr=round(prom/(scatt.R(2)-scatt.R(1)));
klaster_nr(klaster_nr==0)=1;
klaster_nr(klaster_nr>numel(scatt.R))=numel(scatt.R);
rozp=scatt.rel_int(:,I);
Iscat=rozp(klaster_nr);
Iscat(end)=[];
Ix=rozp(klaster_nr(end));
end
clear
clc
load('/home/pracownicy/karpinska/Dokumenty/Ruch_kropli/Dziury_dane/Dane_od_Tiny/Probab_dist_29.mat')
R=size_distribution(10000,Rmin,3,PD);
histogram(R,size(PD,1),'Normalization','probability')
hold on
plot(PD(:,1)*10^(-6),PD(:,2))
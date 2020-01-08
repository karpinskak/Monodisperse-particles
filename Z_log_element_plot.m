%%%%%%%%%% 
% This script calculates real escape time (/divided by tau_p of particle) and
% compares it to the logarithmic element to see where the approximation is
% good
clear
close all
clc
delete(gcp('nocreate'))

%% Declare constants
C = Constants;
%% Load functions
fDIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/Functions/othercolor';
addpath(fDIR)
fun = @(n,m) log((m-sign(n-1))./abs(n-1));

figHandler = figure('Color',[1 1 1]);
width=16;
height=20;
fsize=18;
set(gcf,'Position', [300, 300, 2*560, 2.3*420 ],...
    'paperunits','centimeters',...
    'papersize',[width,height],...
    'InvertHardCopy','off')

% syms x y
% funkcja=@(x) log(abs(x))^2
% int(funkcja,x)
m=[1,10:10:50];
for j=1:numel(m)
    q(j)=integral(@(n) fun(n,m(j)),-m(j),m(j));
    od(j)=integral(@(n) fun(n,m(j)).^2,-m(j),m(j));
    n=linspace(-m(j),m(j),3000);
    for k=1:numel(n)
        pojed(j,k)=fun(n(k),m(j)); 

    end
    plot(n,pojed(j,:),'Linewidth',2,'Color',[0.5 j/numel(m) 0])
    hold on
end

xlim([-m(end),m(end)])
ylim([0,10])
ticzki=sort([-50:10:50,-1]);
ticzki(ticzki==0)=[];
xticks(ticzki)
xlabel('z^\ast_0')
ylabel('log L(Z^\ast,z^\ast_0)')
grid on
h=legend(sprintfc('%d',m));
title(h,'Z^\ast');
set(gca,'FontSize',fsize)

% sr=q./(2*m);
% plot(m,sr)


%% 
% figure(1)
%  m=2:5:1000;
%  p=[10,50,100,150,200,250,300,500,1000,10000,100000];
%  for r=1:numel(p)
%  for j=1:numel(m)
%      n=linspace(-m(j),m(j),p(r));
%      if numel(nonzeros(n==1))~=0
%          srednia(j)=NaN;
%      else
%      czyn(j,:)=log((m(j)-sign(n-1))/(abs(n-1)));
%      srednia(j)=mean( czyn(j,:));
%      end
%  end
%  plot(m,srednia)
%  str{r}=char(num2str(p(r)));
%  hold on
%  end
%  legend(str)
%  
%  
%  function czaskon=wylicz_czas_wylotu(z_0,nu,delta,tau_p,Z,g,teta)
% if Z<=0
%     error("This function takes only Z>0")
% end
% 
% k=(sqrt(1+(nu*tau_p/(delta^2)))).^(-1);
% z_b=g*delta^2*tau_p*cos(teta)/nu;
% tau_2=2*tau_p/(k^(-1)-1);
% czaskon=(log((Z-z_b*sign(z_0-z_b))./abs(z_0-z_b))-log((1+k)/2))*tau_2;
% end
%  
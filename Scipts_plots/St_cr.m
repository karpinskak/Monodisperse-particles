clear
clc
clf

figure(1)
%rr=0:5*10^(-3):1.585;
r=[(0.909542834978651^2+0.111974184899316^2)^(1/2),(0.892269830924517^2+0.154722616396216^2)^(1/2),(0.881782985548095^2+0.175111112665929^2)^(1/2)];
%AA=10^(-6):10^(-4):0.040001;
A=[.008,.01132,.013];
%[r,A]=meshgrid(rr,AA);
%pcolor(fun_fi(r))
Stcr=A./abs(fun_fi(r));
h=surf(A,r,log(Stcr));
set(h, 'EdgeColor', 'none');
colorbar

function wart=fun_fi(r)
part1=(1-exp(-r.^2/2))./r.^2;
part3=exp(-r.^2/2);
wart=part1.*(part1-part3)/(2*pi)^2;
end
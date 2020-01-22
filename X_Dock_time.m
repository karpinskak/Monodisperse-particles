% This script calculates docking time of particles
% Docking time is a time at which particle positioned at r=0 can
% reach r=r_s position or positioned at r=r_s can reach r=0. There is no
% gravity, so there is axial symmetry. Particles start at fi=0. Final
% simulation time is nondimensional. We analyze
% dependence on A/St and A.
 
clear
close all
clc
delete(gcp('nocreate'))
 
tic
%% Declare constants and parameters
Const = Constants;
X0_Parameters
 
 
%% Initial conditions generation
%polozenie = 0, rs, plus epsilon
 
part.init=[];
part.traj=[];
part.traj.X=[];
part.traj.t=[];
 
part.par=[];
part.par.St=[];
part.par.A=[];
part.init.r0_bez=[];
part.init.rfin_bez=[];
 
fi0_bez=0;
r0_bez=zeros(1,numel(A));
rfin_bez=zeros(1,numel(A));
texit=zeros(1,numel(A));
texit2=zeros(1,numel(A));
vr_0=zeros(1,numel(A));
vfi_0=zeros(1,numel(A));
flag=zeros(1,numel(A)); % flag=0 - point docking, flag=1 - in orbit docking
iterations=zeros(1,numel(A));
calc_time=zeros(1,numel(A));
 
%% calc
parpool('local',poolnr)
 
parfor p=1:numel(A)
       
    disp(char(num2str(p),' is on.'))
    tic
    t0=t00;
    tfin=St(p)*nt;
    calc_time_temp=0;
    
    if St(p)/A(p)<Const.St_A_cr
        r0_bez(p)=Const.rs-eps;
        rfin_bez(p)=eps;
        [vr_0(p),vfi_0(p),~]=velocity_field(A(p),Const.rs-eps,0);
    elseif St(p)/A(p)>=Const.St_A_cr
        r0_bez(p)=eps;
        rfin_bez(p)=row_orb_bur(A(p)/St(p))-eps;
        [vr_0(p),vfi_0(p),~]=velocity_field(A(p),eps,0);
        flag(p)=1;
    end
    
    x00=[r0_bez(p),...
        fi0_bez,...
        vr_0(p),...
        vfi_0(p)];
    
    rfin_temp=rfin_bez(p);
    tspan=linspace(t0,tfin,nt);
    
    opt = odeset('Events',@(t, x) EventFunction(t, x, rfin_temp));
    %
    [t1,x1,tev, ~, iev]=ode45(@(t,x)traj_pionowy(t,x,St(p), A(p)),tspan,x00,opt);
    iteracja=0;
    trjx=x1;
    trjt=t1;
        
    while isempty(iev) && calc_time_temp<lim_time % poki nie znaleziono wydarzenia i czas obliczen nie przekracza limitu
        % sa dwie mozliwe przyczyny, dlaczego event sie nie wydarzyl
        
        if (x1(2,1)-x1(1,1))/(rfin_bez(p)-r0_bez(p))>1
            % gdy dt bylo za duze i event zdazyl sie miedzy t(1) a
            % t(2), a program tego nie wychwycil
            tfin=t1(2);
            tspan=linspace(0,tfin,nt);
            [t1,x1,tev, ~, iev]=ode45(@(t,x)traj_pionowy(t,x,St(p), A(p)),tspan,x00,opt);
            trjx=x1;
            trjt=t1;
        elseif (x1(2,1)-x1(1,1))/(rfin_bez(p)-r0_bez(p))<=1
            %dt nie bylo za duze, czyli czas koncowy musial byc za krotki zeby dotrzec do eventu
            t0=tfin;
            pr=(rfin_bez(p)-x1(end,1))/mean(x1((end-10):end,3));
            if pr<0 || pr>5*tfin
                tfin=4*tfin;
            else
            tfin=t0+1.1*pr;
            end
            tspan=linspace(t0,tfin,nt);
            if ~isnan(tfin)
                x0=x1(end,:);
                [t1,x1,tev, ~, iev]=ode45(@(t,x)traj_pionowy(t,x,St(p),A(p)),tspan,x0,opt);
            else
                tev=NaN;
            end
            trjx=[trjx;x1];
            trjt=[trjt;t1];
        end
        iteracja=iteracja+1;
        calc_time_temp=toc;
    end
    
    if isempty(iev) && calc_time_temp>=lim_time
        texit(p)=1i*tfin;
    elseif ~isempty(iev)
        texit(p)=tev;
        % now calculate it with enough detail, take small timestep adjusted to particle size
        % it gives the range around tev that will become tspan for detailed traj calc
        nr=2;
        tdif=trjt(end)-trjt(end-nr);
        while tdif<0.5*St(p)
            nr=nr+1;
            tdif=trjt(end)-trjt(end-nr);
        end
        tspanev=trjt(end-nr):0.2*St(p):(trjt(end)+tdif);
        x0ev=trjx(end-nr,:);
        [t1,x1,tev2, ~, iev2]=ode45(@(t,x)traj_pionowy(t,x,St(p),A(p)),tspanev,x0ev,opt);
        if ~isempty(tev2)
            texit2(p)=tev2;
        end
    end
    
    iterations(p)=iteracja;
    calc_time(p)=calc_time_temp;
    trX{p}=trjx;
    trt{p}=trjt;
    disp(char('Calculation for p=',num2str(p),' done.'))
    
    if sign(rfin_bez(p)-r0_bez(p))~=sign(trjx(end,1)-trjx(1,1))
        disp(char('Sth is wrong with particle direction of motion for p=',num2str(p),'.'))
        texit(p)=NaN;
    end
    % the end of calculation for p
end
 
delete(gcp('nocreate'))
%plot(reshape(A,size(Aa)),reshape(St,size(Aa)),reshape(texit,size(Aa)));
 
%% saving
for p=1:numel(A)
    part(p).par.St=St(p);
    part(p).par.A=A(p);
    part(p).init.r0_bez=r0_bez(p);
     part(p).init.rfin_bez=rfin_bez(p);
     part(p).init.fi0_bez=fi0_bez;
     part(p).traj.X=cell2mat(trX(p));
     part(p).traj.t=cell2mat(trt(p));
end
texitt=texit; %% changed 22.01.2020
texitt(texitt==0)=NaN;
toc
save([loadDIR name])
 
function dx = traj_pionowy(t,x,St,A)
dx=[x(3);...
    x(4);...
    -A*x(1)/St+x(1)*(x(4))^2-x(3)/St;...
   (1-exp(-x(1)^2/2))/(2*pi*St*x(1)^2)-2*x(3)*x(4)/x(1)-x(4)/St];
end
 
function [value,isterminal,direction] = EventFunction(t,x,P)
value=x(1)-P; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end
 
function [ r0 ] = row_orb_bur( par )
%ROWNANIE_ORBITY wyznacza promieĹ„ orbity stacjonarnej wiru Burgersa bez
%wpĹ‚ywu grawitacji biorÄ…c parametr m=A/St.
 
syms r
r0_sym=solve((par)^(1/2)*r^2-(1-exp(-r^2/2))/(2*pi)==0);
 
r0=abs(double(r0_sym(2)));
 end
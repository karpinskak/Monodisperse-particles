

if par_set==0
    P=part(1).par.gamma*part(1).par.A*part(1).par.delta^2/Const.nu;
    delta_t=0.1*part(1).par.tau_f;
else
    M1=part(1).par.gamma*part(1).par.delta^4*part(1).par.A^2./(tau_p*(Const.nu)^2);
    M2=part(1).par.delta^2*part(1).par.A./(Const.nu*tau_p);
    M3=Const.g*part(1).par.delta^3*sin(part(1).par.teta)*part(1).par.A^2/((Const.nu)^2);
    M4=Const.g*part(1).par.delta^3*cos(part(1).par.teta)*part(1).par.A^2/((Const.nu)^2);
    delta_t=0.2*part(1).par.tau_p;
   
end


tglobal=tstart:delta_t:(T+tstart);

s =  isempty(gcp('nocreate'));
 if s==1
 parpool('local',poolnr)
 end
Indeks=zeros(l_krop,2);

parfor l=1:l_krop
    
    %
    M1r=M1(l);
    M2r=M2(l);
    %
    [~,I0]=min(abs(part(l).init.t0-tglobal));  % punkty czasu globalnego najbliższe czasu początkowego i końcowego
    [~,Ik]=min(abs(part(l).init.tfin-tglobal));

    id=[0,0];
    if part(l).init.t0-tglobal(I0)>10^(-8)  % liczba ta wyraza najmniejsza mozliwa dopuszczalna roznice miedzy punktami ruchu
        I0=I0+1;
    elseif abs(part(l).init.t0-tglobal(I0))<10^(-8)
        id(1)=1;
    end
    
    if part(l).init.tfin-tglobal(Ik)<-10^(-8)
         Ik=Ik-1;
    elseif abs(part(l).init.tfin-tglobal(Ik))<10^(-8)
        id(2)=1;
    end
    
    if I0<=Ik
        
        tzakres=tglobal(I0:Ik);
        if id(1)==0
            tzakres=[part(l).init.t0,tzakres];
        end
        if id(2)==0
            tzakres=[tzakres,part(l).init.tfin];
        end
        tspan=tzakres/part(1).par.tau_f;
        
        
        % rescale to nondimensional
        x0=[part(l).init.r0/part(1).par.delta,...
            part(l).init.fi0];
        z0=part(l).init.z0/part(1).par.delta;
        
        %trajectory calculation
        if par_set~=0
             x0=[x0,...
                (-part(1).par.gamma*part(l).init.r0/2)/part(1).par.v_f,...
                0];
            z0=[z0,0];
            [t1,x1]=ode45(@(t,x)trajektoria_polarne(t,x,M1r,M2r,M3),tspan,x0);
            [~,z1]=ode45(@(t,z)trajektoria_z(t,z,M1r,M2r,M4),tspan,z0);
        else
            [t1,x1]=ode45(@(t,x)traj_pol_skalar(t,x,P),tspan,x0);
            [~,z1]=ode45(@(t,z)traj_z_skalar(t,z,P),tspan,z0);
        end
        
        % rescale to real dimensions
        
        r=x1(:,1)*part(1).par.delta;
        fi=x1(:,2);
        z=z1(:,1)*part(1).par.delta;
   
      
       first=2-id(1);
       last_normal=numel(tzakres)-1+id(2);
       last_nonzero=find(r<=part(1).par.D,1,'last');
       if isempty(last_nonzero)~=1
       last=min(last_normal,last_nonzero);
       else
           last=last_normal;
       end
       
        t=t1(first:last)*part(1).par.tau_f;
        
        %Przejscie do wspolrzednych kartezjanskich na potrzeby wyrysowania wykresu
        X=r(first:last).*cos(fi(first:last));
        Y=r(first:last).*sin(fi(first:last));
        Z=z(first:last);       

        trt{l}=t;
        trX{l}=X;
        trY{l}=Y;
        trZ{l}=Z;
        
        Iknew=find(abs(tglobal-t(end))<10^(-8));
        Indeks(l,:)=[I0,Iknew];
        
    else
        trt{l}=[];
        trX{l}=[];
        trY{l}=[];
        trZ{l}=[];
        Indeks(l,:)=[0,0];     
    end
end
toc
disp('Trajectory calculation itself finished.')

for  l=1:l_krop
    part(l).traj.t=cell2mat(trt(l));
    part(l).traj.X=cell2mat(trX(l));
    part(l).traj.Y=cell2mat(trY(l));
    part(l).traj.Z=cell2mat(trZ(l));
end

clear trt trX trY trZ l
if par_set~=0
    clear M1 M2 M3 M4
else
    clear P
end
save([loadDIR '/Trajectories.mat'],'part','delta_t','tglobal','l_krop','-v7.3')
toc
disp('Trajectory data saved.')

 parfor k=1:numel(tglobal)
    numerki=[];
    for l=1:l_krop
    if k>=Indeks(l,1)&& k<=Indeks(l,2)
        numerki=[numerki,l];
    end
    end
        drop_in_time{k}=numerki;
 end   

save([loadDIR '/Trajectories.mat'],'drop_in_time','-append');


function dz=trajektoria_z(t,y,M1r,M2r,M4r)
dz=[y(2); M1r*y(1)-M2r*y(2)-M4r];
end

function dx=trajektoria_polarne(t,x,M1r,M2r,M3r)
dx=[x(3);x(4); -M1r*x(1)/2-M2r*x(3)-M3r*sin(x(2))+x(1)*(x(4))^2; M2r*(1-exp(-x(1)^2/2))/(2*pi*x(1)^2)-M3r*cos(x(2))/x(1)-2*x(3)*x(4)/x(1)-M2r*x(4)];
% w tym wzorze jest bład: zamiast + jest * przed ostatnim członen
% pierwszego rownania
end

function dx=traj_pol_skalar(t,x,P)
dx=[-P*x(1)/2;...
    (1-exp(-x(1)^2/2))/(2*pi*x(1)^2)];
end

function dz=traj_z_skalar(t,z,P)
dz=P*z;
end

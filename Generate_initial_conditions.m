switch type
    case 0
        switch new_sim
            case 0
                [polozenia_0]=positions_t0(max_R,part.par.D,part.par.Z,n,tstart);
                [polozenia_1]=positions_rand(max_R,part.par.D,part.par.Z,n,part.par.gamma,tstart,T,poolnr);
                polozenia=[polozenia_0;polozenia_1];
                clear polozenia_0 polozenia_1
                
                %% Final time of particle in the vortex
                switch dispersity
                    case 0
                        Rall = size_distribution(size(polozenia,1),part.par.Rmin,dispersity,part.par.R);
                    case 1
                        Rall = size_distribution(size(polozenia,1),part.par.Rmin,dispersity,part.par.Rmax);
                    case 2
                        Rall = size_distribution(size(polozenia,1),part.par.Rmin,dispersity,part.par.R,part.par.Rdev);
                    case 3
                        Rall = size_distribution(size(polozenia,1),part.par.Rmin,dispersity,part.par.ProbDist);
                end
                [part_par,tau_p]=nowe_parametry_kropel(Const,Rall,part.par.delta,part.par.A,part.par.teta,part.par.gamma,part.par.tau_f);
                czaskon=wylicz_czas_wylotu(polozenia(:,4),part.par.tau_f,part.par.gamma,part.par.Z,part_par.tau_p,part_par.z_b);
                polozenia(:,5)=((czaskon+polozenia(:,1))<=T+tstart).*(czaskon+polozenia(:,1))+((czaskon+polozenia(:,1))>T+tstart).*(T+tstart);
                clear czaskon
                
            case 1
                load([loadPosDIR,'/Init.mat'],'polozenia','part_par') % load file with starting positions
                % check if this is for the same D,Z,n,max_R, gamma, T
                
                %     case 2
                %         load([loadDIR,'coscos.mat']) % load last positions of the particles
                %         % (at tstart) as new starting positions
                %         % check if this is for the same D,Z,n,max_R, gamma
                %         % polozenia_00=()
                %         [polozenia_2,l_krop]=positions(max_R,part.par.D,part.par.Z,n,part.par.gamma,tstart,T);
                %         polozenia=[polozenia_00;polozenia_2];
                %         % polacz pozniej polozenia_1 i 2
        end
    case 1
        polozenia=positions_even(part.par.D,n,tstart,T);
        switch dispersity
            case 0
                Rall = size_distribution(size(polozenia,1),Rmin,dispersity,part.par.R);
            case 1
                Rall = size_distribution(size(polozenia,1),Rmin,dispersity,part.par.R,part.par.Rmax);
            case 2
                Rall = size_distribution(size(polozenia,1),Rmin,dispersity,part.par.R,part.par.Rdev);
            case 3
                Rall = size_distribution(size(polozenia,1),Rmin,dispersity,part.par.ProbDist);
        end
        [part_par,tau_p]=nowe_parametry_kropel(Const,Rall,part.par.delta,part.par.A,part.par.teta,part.par.gamma,part.par.tau_f);
end

l_krop=numel(polozenia(:,1));
%% save
part.init=[];
part.init.t0=[];
part.init.r0=[];
part.init.fi0=[];
part.init.z0=[];
part.init.tfin=[];
for p=1:l_krop
    part(p).init.t0=polozenia(p,1);
    part(p).init.r0=polozenia(p,2);
    part(p).init.fi0=polozenia(p,3);
    part(p).init.z0=polozenia(p,4);
    part(p).init.tfin=polozenia(p,5);
    part(p).new_par.tau_p=part_par.tau_p(p);
    part(p).new_par.tau_g=part_par.tau_g(p);
    part(p).new_par.St=part_par.St(p);
    part(p).new_par.Sv=part_par.Sv(p);
    part(p).new_par.Fr=part_par.Fr(p);
    part(p).new_par.B=part_par.Fr(p);
    part(p).new_par.z_b=part_par.z_b(p);
    part(p).new_par.R=Rall(p);
end
clear p
part(1).par.l_krop=l_krop;
if new_sim==0
    save([loadPosDIR '/Init.mat'],'polozenia','part_par','tau_p')
end

%% functions
function [polozenia]=positions_even(D,n,tstart,T)
ile_kaw=floor(2*D*(n^(1/3)-1))+mod(ceil(2*D*(n^(1/3)-1)),2);
Xlin=-D:2*D/ile_kaw:D;
Ylin=-D:2*D/ile_kaw:D;
[X,Y]=meshgrid(Xlin,Ylin);
R=(X.^2+Y.^2).^(1/2);
FI=atan2(Y,X);
pol_r=reshape(R,[numel(Xlin)*numel(Ylin),1]);
pol_fi=reshape(FI,[numel(Xlin)*numel(Ylin),1]);
if numel(pol_r)~=numel(pol_fi)
    error('Cos nie tak z losowaniem jednorodnym')
end
polozenia=zeros(numel(pol_r),5);
polozenia(:,2)=pol_r;
polozenia(:,3)=pol_fi;
polozenia(:,4)=0;
polozenia(:,5)=tstart+T;
polozenia(:,1)=tstart;

end

function [polozenia]=positions_rand(max_R,D,Z,n,gamma,tstart,T,poolnr)

% generate grid
delta_r=2*max_R;

u=gamma*D/2;
delta_fi=delta_r/D;
V=delta_fi*D*delta_r^2; % objetosc elementarnego szescianiku siatki
Nv=n*V; % liczba kropel w objetosci elementarnej (prawdopodobienstwo)

q=floor(2*pi/delta_fi)+1; % liczba szescianikow
p=floor(2*Z/delta_r)+1;   %

tau=delta_r/u; % czas, w jakim kropla opusci elementarna komorke na powierzchnii

% declare approx. nr of 
l_kr_los= 5*ceil(p*q*Nv);% maksymalna przewidywana l wylosowanych w jednym losowaniu kropel
t0_wek=tstart+tau:tau:tstart+T;
l_los=numel(t0_wek);
pol_fi_m=zeros(l_kr_los,l_los); % w wersji macierzy a nie kolumny
pol_z_m=zeros(l_kr_los,l_los);
l_kropel=zeros(1,l_los);
czasy=zeros(l_kr_los,l_los);

s =  isempty(gcp('nocreate'));
if s==1
parpool('local',poolnr)
end

parfor k=1:numel(t0_wek)
   
M=rand(p,q);
krople=double(M<Nv);
[row,col]=find(krople);
l_krop=numel(row);
l_kropel(k)=l_krop;
t0=t0_wek(k);
czasik=zeros(l_krop,1)+t0;
if l_krop~=0
czasy(:,k)=[czasik;zeros((l_kr_los-l_krop),1)];
pol_fi=single((col-1)*delta_fi);
pol_fi_m(:,k)=[pol_fi;zeros((l_kr_los-l_krop),1)];
pol_z=single(-Z+delta_r*(row-1));
pol_z_m(:,k)=[pol_z;zeros((l_kr_los-l_krop),1)];
end
end

calk_l_krop=sum(l_kropel);
polozenia=zeros(calk_l_krop,4);
dot_l_krop=0;
for k=1:numel(t0_wek)
  
polozenia((dot_l_krop+1):(dot_l_krop+l_kropel(k)),1)=czasy(1:l_kropel(k),k);
polozenia((dot_l_krop+1):(dot_l_krop+l_kropel(k)),3)=pol_fi_m(1:l_kropel(k),k);
polozenia((dot_l_krop+1):(dot_l_krop+l_kropel(k)),4)=pol_z_m(1:l_kropel(k),k);

 dot_l_krop=l_kropel(k)+dot_l_krop;% dotychczasowa liczba kropel
end
polozenia(:,2)=D-0.5*delta_r;

l_krop=numel(polozenia(:,1));
l_krop_teoret=l_los*n*pi*(D^2-(D-delta_r)^2)*2*Z;

if (l_krop_teoret-l_krop)/l_krop_teoret>0.05
    error('Wrong nr of particles generated')
end
end

function [polozenia_0]=positions_t0(max_R,D,Z,n,tstart)
%% Find initial positions of particles for t=0

delta_xyz=6*max_R;
Nv=n*delta_xyz^3; % nr of particles in the single cube = probability
% first generate them in the whole cuboid
p=floor(2*D/delta_xyz)+1; % nr of cubes along X and Y
s=floor(2*Z/delta_xyz)+1;  % nr of cubes along Z

M=rand(p,p,s);
krople=double(M<Nv);
ind=find(krople);
[row,col,page]=ind2sub(size(krople), ind);


% put particle in the middle of the cube
pol_x_all=-D+(row-1)*delta_xyz+0.5*delta_xyz;
pol_y_all=-D+(col-1)*delta_xyz+0.5*delta_xyz;
pol_z_all=-Z+delta_xyz*(page-1)+0.5*delta_xyz;

% cut out particles outside the cyllinder
pol_r_all=(pol_x_all.^2+pol_y_all.^2).^(0.5);
pol_fi_all=atan2(pol_y_all,pol_x_all);

pol_fi=pol_fi_all.*(pol_r_all<=D); 
pol_r=pol_r_all.*(pol_r_all<=D);
pol_z=pol_z_all.*(pol_r_all<=D);

%% find artificial particles at (0,0,0)
polozenia_0=zeros(numel(nonzeros(pol_r)),4);

war=pol_fi.*pol_r.*pol_z;
war2=abs(pol_r)+abs(pol_fi)+abs(pol_z);

if numel(nonzeros((war==0).*(war2~=0)))~=0
    disp('Wylosowano element z jednym zerem w polozeniu')
    [zerowy_el,zer_el]=find((war==0).*(war2~=0));
    error('Blad')
else
    polozenia_0(:,2)=nonzeros(pol_r);
    polozenia_0(:,3)=nonzeros(pol_fi);
    polozenia_0(:,4)=nonzeros(pol_z);
    disp('There are no elements with zero in the position')
end

polozenia_0(:,1)=tstart;
l_krop=numel(polozenia_0(:,1));
l_krop_teoret=n*pi*D^2*2*Z;
if (l_krop_teoret-l_krop)/l_krop_teoret>0.05
    error('Wrong nr of particles generated - t0')
end
end

function [new_par,tau_p]=nowe_parametry_kropel(C,Rall,delta,A,teta,gamma,tau_f)

tau_p=(2*C.ro_p*Rall.^2)/(9*C.nu*C.ro_a);
new_par.R=Rall;
new_par.tau_p=tau_p;
new_par.St=new_par.tau_p*C.nu/(delta^2*A);
new_par.Sv=C.g*delta*A*new_par.tau_p*sin(teta)/(C.nu);
new_par.tau_g=tau_f./new_par.Sv;
new_par.Fr=new_par.St./new_par.Sv;
new_par.B=(C.rs/(2^8*pi^3))*(new_par.St*A^2./new_par.Sv);
new_par.z_b=C.g*new_par.tau_p*cos(teta)/gamma;
end

function czaskon=wylicz_czas_wylotu(polozenia_z,tau_f,gamma,tau_p,Z,z_b)

lambda1=0.5*tau_f.*(-1+(4*gamma.*tau_p+1).^(0.5))./tau_p;
lambda2=0.5*tau_f.*(-1-(4*gamma.*tau_p+1).^(0.5))./tau_p;
czaskon=abs(tau_f.*log((lambda2-lambda1).*(Z.*sign(polozenia_z-z_b)-z_b)./((polozenia_z-z_b).*lambda2))./lambda1);
end
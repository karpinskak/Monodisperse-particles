clear
close all
clc
delete(gcp('nocreate'))

DIR='/home/pracownicy/karpinska/Dokumenty/Praca_doktorska_analizy/Monodisperse_particles/';

%% Load constants and functions
addpath(DIR)
Const=Constants;
fDIR=[DIR '/Functions'];
addpath(fDIR)
%
nproc=3;

%% choose input data type
dane=2; % 1 - if you have meshgrid of St and A parameters, 2- list of St/A,
        % 3 - meshgrid of R and A for given delta
        % 4 - meshgrid of delta and A
if dane==2 % choose
    opcja=3; % 1 - plot r0 vs St/A, 2- plot r0 vs A/St 3- inne
end

% plot parameters
fsize=16;
width=16;
height=20;


%% Input data and calculation of param=St/A;

switch dane
    case 1  %give matrix of St and A
        ST=[0.0001:0.0001:0.001,0.0011:0.001:0.101,0.111:0.01:1];
        AA=[0.0001:0.0001:0.0081];
        [A,St]=meshgrid(AA,ST);
        param=St./A;
        r0=zeros(size(A));
        
    case 2 % give parameter values list
        param=[(Const.St_A_cr+10^(-20)):1:(10+Const.St_A_cr), (10+Const.St_A_cr):10:10000];
        r0=zeros(size(param));
        
    case 3 % give matrix of R and A for given delta
        
        delta=0.005; %[m]
        skala=100;
        RR=(1:0.5:20)*10^(-6); %[m]
        AA=0.0001:0.00005:0.005; %[1]
        [R,A]=meshgrid(RR,AA);
        St=zeros(size(A));
        
        for k=1:numel(RR)
            for l=1:numel(AA)
                [par]=wylicz_param(Const,2,R(l,k),delta,0,A(l,k),0,0);
                St(l,k)=par.St;
            end
        end
        clear par
        param=St./A;
        r0=zeros(size(A));
end

%% calculate orbits
parpool('local',nproc)
tic
switch dane
    
    case 1
        for j=1:numel(AA)
            r0_temp=zeros(size(ST))';
            parfor k=1:numel(ST)
                ptemp=param(k,j);
                if ptemp>Const.St_A_cr
                    r0_temp(k)=Row_orb_bur(ptemp);
                    
                else
                    r0_temp(k)=NaN;
                end
            end
            r0(:,j)=r0_temp;
        end
        toc
        
    case 2 
        parfor k=1:numel(param)
            ptemp=param(k);
            disp(['Calculating orbit for k=',num2str(k)])
            if ptemp>Const.St_A_cr
                r0(k)=Row_orb_Bur(ptemp);
            else
                r0(k)=NaN;
            end
            disp(['Done for k=',num2str(k)])
        end
        toc
        
    case 3
        for j=1:size(param,1)
            r0_temp=zeros(1,size(param,2))';
            parfor k=1:size(param,2)
                ptemp=param(j,k);
                if ptemp>Const.St_A_cr
                    r0_temp(k)=Row_orb_Bur(ptemp);
                    
                else
                    r0_temp(k)=NaN;
                end
            end
            r0(j,:)=r0_temp;
        end
        toc
end

%% plot

switch dane
    
    case 1
        surf(A,St,r0)
        colorbar
        xlabel('A')
        ylabel('St')
        xlim([0.0001 1])
        ylim([0.0001 0.03])
        
    case 2
        close all
        figure
        set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
            'paperunits','centimeters',...
            'papersize',[width,height],...
            'InvertHardCopy','off')
        
        switch opcja
            case 2
                % yyaxis left
                line(param.^(-1/2),r0,'Color','r')
                ax1 = gca; % current axes
                ax1.XColor = 'r';
                ax1.YColor = 'r';
                ticzki=get(gca, 'XTick');
                ticzki(ticzki==0)=[];
                ticz=unique([1/Const.St_A_cr, ticzki]);
                set(gca, 'XTick', ticz,...
                    'FontSize',fsize);
                xlabel('\sqrt{A/St}')
                ylabel('r^+_{orb}')
                hold on
                ax1_pos = ax1.Position; % position of first axes
                ax2 = axes('Position',ax1_pos,...
                    'XAxisLocation','top',...
                    'YAxisLocation','right',...
                    'FontSize',fsize,...
                    'Color','none')
                %'Xlim',[5 9.21],...
                
                % yyaxis right
                %omega=sqrt(param.^(-1));
                %line(param,omega.*r0.^2,'Parent',ax1,'LineWidth',2)
                %ylabel('\omega_{orb}')
                %hold on
                
                line(log(param.^(-1/2)),r0,'Parent',ax2,'Color','k')
                xlabel('log \sqrt{A/St}')
                ylabel('r^+_{orb}')
                
            case 1
                yyaxis left
                line(param.^(0.5),r0,'Color','k','LineWidth',2)
                ax1 = gca; % current axes
                ax1.XColor = 'k';
                ax1.YColor = 'k';
                ticzki=get(gca, 'XTick');
                ticzki(ticzki==0)=[];
                ticz=unique([floor(Const.St_A_cr), ticzki]);
                set(gca, 'XTick', ticz,...
                    'FontSize',fsize);
                xlabel('$\sqrt{St/A}$','interpreter','latex')
                ylabel('$r^+_{orb}$','interpreter','latex')
                hold on
                
                yyaxis right
                omega=sqrt(param.^(-1));
                line(param.^(0.5),omega.*r0.^2,'Parent',ax1,'LineWidth',2)
                ylabel('$\omega_{orb}$','interpreter','latex')
                
                ax1_pos = ax1.Position; % position of first axes
                ax2 = axes('Position',ax1_pos,...
                    'XAxisLocation','top',...
                    'YAxisLocation','left',...
                    'Xscale','log',...
                    'FontSize',fsize,...
                    'Color','none');
                %'Xlim',[1.48*10^(2),10^(4)],...
                
                line(param.^(0.5),r0,'Parent',ax2,'Color','b','LineWidth',2)
                set(gca,'XColor','b')
                xlabel('$\sqrt{St/A}$','interpreter','latex')
                ylabel('$r^+_{orb}$','interpreter','latex')
                box on
                hold on
                plot(4*pi+zeros(size(r0)),r0,'-.','Color',[0.65,0.65,0.65],'LineWidth',1)
                
            case 3
                yyaxis left
                line(param.^(0.5),r0.^2,'LineWidth',2)
                %loglog(param.^(0.5),r0,'LineWidth',2)
                ylabel('$r^{+ 2}_{orb}$','interpreter','latex')
                hold on
                line(param.^(0.5), param.^(1/2)/(2*pi),'LineStyle','--','LineWidth',1)
                line(param.^(0.5), 4-16*pi*param.^(-0.5),'LineStyle','--','LineWidth',1)
                
                yyaxis right
               omega=sqrt(param.^(-1));
                %loglog(param.^(0.5),omega.*r0.^2,'LineWidth',2)
                line(param.^(0.5),omega.*r0.^2,'LineWidth',2)
                ylabel('$\omega_{orb} r^{+ 2}_{orb}$','interpreter','latex')
                
                xlabel('$\sqrt{St/A}$','interpreter','latex')
                line(param.^(0.5),1/(2*pi)+zeros(size(param)),'LineStyle','--','LineWidth',1)
                set(gca,'FontSize',fsize)
        end
        
        
        
    case 3
        close all
        figure
        set(gcf,'Position', [640, 300, 2*560, 2*420 ],...
            'paperunits','centimeters',...
            'papersize',[width,height],...
            'InvertHardCopy','off')
        
        pcolor(A,R*10^6,r0*delta*skala)
        hold on
        
        St_cr=16*pi^2*AA;
        tau_p_cr=delta^2*AA.*St_cr*Const.nu^(-1);
        R_cr=(9*Const.nu*tau_p_cr*Const.ro_a/(2*Const.ro_p)).^(1/2);
        clear St_cr tau_p_cr
        plot(AA,R_cr*10^6,'Color',[0 0 0],'LineWidth',1)
        xlabel('A')
        ylabel('R [um]')
        h = colorbar;
        ylabel(h, 'r_{orb} [cm]')
        shading interp
        ylim([0 20])
        xlim([0.0001 0.005])
        set(gca,'FontSize',fsize)
end

%% saving
delete(gcp('nocreate'))
switch dane
    case 2
        nazwa=char(strcat(DIR, 'Plots/orbit_radius/','orbit_vs_St_A_logSt_A'));
        plik=char(strcat(DIR, 'Results_orbits/Orbits_parameter_vs_St_dev_A.mat'));
    case 1
        nazwa=char(strcat(DIR, 'Plots/orbit_radius/','orbit_surf_St_and_A'));
        plik=char(strcat(DIR, 'Results_orbits/Orbits_parameter_vs_St_and_A.mat'));
    case 3
        nazwa=char(strcat(DIR, 'Plots/orbit_radius/','orbit_vs_R_and_A_delta',strrep(num2str(delta*skala),'.',''),'cm'));
        plik=char(strcat(DIR, 'Results_orbits/Orbits_parameter_vs_R_and_A_delta',strrep(num2str(delta*skala),'.',''),'cm.mat'));
end
saveas(gcf,nazwa,'fig')
saveas(gcf,nazwa,'png')

save(plik)


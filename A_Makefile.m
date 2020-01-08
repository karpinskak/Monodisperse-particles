clear
close all
clc
delete(gcp('nocreate'))
tic

%% Declare constants
Const = Constants;
A0_Parameters

%% Load functions
%fDIR = [DIR 'Trajectories/Functions'];
fDIR=[DIR '/Functions'];
addpath(fDIR)

%% Create directory if not existing
if isfolder(loadDIR)~=1
    mkdir(loadDIR)
end

%% Calculate all parameters needed
[parameters]=wylicz_param(Const, par_set,a,b,c,d,k,l);
clear a b c d k l
part.par=parameters;
switch dispersity
    case 2
        part.par.Rdev=Rdev;
        part.par.ProbDist=PD;
    case 1
        part.par.Rmax=Rmax;
    case 3
        part.par.ProbDist=PD;
end
part.par.Rmin=Rmin;

save([loadDIR '/Parameters.mat'],'max_R','n','T','tstart','parameters','dispersity','par_set','type','new_sim')
clear parameters
fprintf('Parameters generated')

toc
%% Create initial conditions
Generate_initial_conditions
fprintf('Initial conditions generated')

toc
%% Calculate trajectories
Trajectories
fprintf('Trajectories calculated')

toc
delete(gcp('nocreate'))

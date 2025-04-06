%%  Perform DREAM algorithm to CoLM parameter estimation
clc
clear

%% Problem settings defined by user
DREAMPar.d = 6;                          % Dimension of the problem
DREAMPar.lik = 2;                       % Model output is simulation: Gaussian likelihood function

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'fold';         % Explicit boundary handling
Par_info.min = [10,   0.25,  2.5,   0.05,   2.5,   50;];   % If 'latin', min values
Par_info.max = [200,  0.75,  7.5,   0.08,   7.5,   500];   % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'CoLM';
observed_LE=importdata('./Observations/LE_benchmark.txt');   % the filepath of the ynthetic observation of LE
observed_NEE=importdata('./Observations/NEE_benchmark.txt'); % the filepath of the ynthetic observation of NEE
observed_RSM=importdata('./Observations/RSM_benchmark.txt'); % the filepath of the ynthetic observation of RSM
plugin.observed_LE=observed_LE;
plugin.observed_NEE=observed_NEE;
plugin.observed_RSM=observed_RSM;

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = DREAMPar.d;                    % Number of Markov chains
        DREAMPar.T = 10000;                         % Number of generations
    case {'dream_zs','dream_dzs'}
        DREAMPar.N = 3;                             % Number of Markov chains
        DREAMPar.T = 10000;                        % Number of generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.min = -50*ones(1,DREAMPar.d);          % Min value for discrete sampling
    Par_info.max =  50*ones(1,DREAMPar.d);          % Max value for discrete sampling
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % Number of discrete steps
end

%% Optional settings
options.modout = 'no';                % Return model (function) simulations of samples (yes/no)?
options.parallel = 'yes';              % Run each chain on a different core
options.IO='no';
options.save = 'yes';                 % Save workspace DREAM during run

%% Run DREAM package
[chain,output,FX,Z,logL] = DREAM_package(method,Func_name,DREAMPar,Par_info,[],options,plugin);
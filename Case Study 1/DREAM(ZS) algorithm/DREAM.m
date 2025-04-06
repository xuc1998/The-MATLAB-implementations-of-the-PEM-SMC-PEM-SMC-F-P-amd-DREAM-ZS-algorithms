% ---------------------------- Check the following two papers ----------------------------- %
%                                                                                           %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman    %
%       (2009), Accelerating Markov chain Monte Carlo simulation by differential evolution  %
%       with self-adaptive randomized subspace sampling, International Journal of Nonlinear %
%       Sciences and Numerical Simulation, 10(3), 271-288.                                  %
%   Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain with      %
%       snooker updater and fewer chains, Statistics and Computing,                         %
%       10.1007/s11222-008-9104-9.                                                          %
%                                                                                           %
% ----------------------------------------------------------------------------------------- %

%% Problem settings defined by user
DREAMPar.d = 2;                       % Dimension of the problem
DREAMPar.thinning = 1;                 % Only store each 10th sample
DREAMPar.lik = 2;                       % Model output is log-likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.min = [-1 -1]; % If 'latin', min values
Par_info.max = [10 10];  % If 'latin', max values
Par_info.boundhandling = 'none'; 

%% Define name of function (.m file) for posterior exploration
Func_name = 'target';

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
options.modout = 'yes';                % Return model (function) simulations of samples (yes/no)?
options.parallel = 'no';              % Run each chain on a different core
options.save = 'yes';                 % Save workspace DREAM during run


%% Run DREAM package
[chain,output,FX,Z,logL] = DREAM_package(method,Func_name,DREAMPar,Par_info,[],options);
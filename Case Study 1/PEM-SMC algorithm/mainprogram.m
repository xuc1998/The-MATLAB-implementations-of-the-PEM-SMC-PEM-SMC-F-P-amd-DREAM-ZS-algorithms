clc
clear

%% Configuration
Np=3000;                      % number of particles
S=800;                        % number of iterations
bound=[-1 -1            
       10 10];                   % parameter boundary

%% run the PEM-SMC algorithm
[paramter_iteration]=PEM_sampler(Np,S,bound); 
parameter=paramter_iteration(:,:,end);

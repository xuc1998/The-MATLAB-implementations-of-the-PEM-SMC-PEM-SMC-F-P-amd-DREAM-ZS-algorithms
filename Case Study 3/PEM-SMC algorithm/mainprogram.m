clc
clear 
%% Configurations 

Np=200;                      % number of particles
S=100;                      % number of iterations

%% estimating the CoLM model parameters using observed datasets
%  parameter ranges: vmax25 porsl_s bsw_s effcon bsw_d phi0_d  
bound=[10,  0.25,  2.5,   0.05,    2.5,   50; 
       200, 0.75,  7.5,   0.08,    7.5,   500];  
[parameter_iteration,Neff]=PEM_sampler(Np,S,bound); 



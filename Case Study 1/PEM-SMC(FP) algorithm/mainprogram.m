clc
clear

%% Configuration
Np=3000;
S=800;
bound=[-1 -1
       10 10];
[parameter_iteration]=PEM_SMC_FP(Np,S,bound);

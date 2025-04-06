function [parameter_iteration] = PEM_SMC_FP(Np,S,bound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PEM-SMC (F&P) algorithm
% Data: Nov 13, 2024
% Authors: Xu Cong
% Location: Lanzhou University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1. internal settings
dem=size(bound,2);                  % the dimesion of the parameters 
parameter_iteration=nan(Np,dem,S);  % matric used to save samplers duirng the iteration


%% 2. Initialize
%  2.1 Initial population and weight
for i=1:Np
    for d=1:dem
        parameter_old(i,d)=bound(1,d)+(bound(2,d)-bound(1,d))*rand;
    end 
    w(i,1)=1/Np;
end

% 2.2 Prior distribution:  Drawing points uniformly from the prior range
L_P0=0;
LB=bound(1,:);
UB=bound(2,:); 
for i=1:length(LB)
    L_P0=L_P0-log(UB(i)-LB(i));
end

% 2.3 generation of the exponential sequence {beta_s}
m=7e-16; 
x=log(1/m)/log(S);
[betas]=sequencesGen(S,m,x);
betas(S+1:S+5)=1;        % the last five values are forced to be 1

% 2.4 Parallel settings
parpool('local',20); % open parallel pool

%% 3. SMC evolution
parameter_iteration(:,:,1)=parameter_old;
for stage=2:S
    % 3.1 Reweighting
    betas_new=betas(stage);
    betas_old=betas(stage-1);   
    for k=1:Np
        targtet1=target(parameter_old(k,:)); 
        alp=betas_new* targtet1;     % log  of target distribution
        prix=(1-betas_new)*L_P0;                        % log  of prior distribution    
        
        alp_old=betas_old*targtet1;
        prix_old=(1-betas_old)*L_P0;           
        cal_weight(k,1)=(alp+prix)-(alp_old+prix_old);  
    end
    w=w.*exp(cal_weight);
    w=w./(sum(w));% Normalize

    % 3.2 Resampling
    ind=ResampSys(w, Np);
    parameter_new=parameter_old(ind,:);    
    w=repmat((1/Np),Np,1);  

    % 3.3 Random walk
    parfor i=1:Np
        parameter_random=randw(parameter_new(i,:),LB,UB);
        alp_old=betas_new*target(parameter_new(i,:));
        alp_media=betas_new*target(parameter_random);    
         % M-H accept
        ratio=min([1,exp(alp_media-alp_old)]);
        u=rand;
        if u<ratio
            parameter_new(i,:)=parameter_random;  
        else
            parameter_new(i,:)=parameter_new(i,:); 
        end        
    end    

    % 3.4 Mutation
    parameter_new_updated=zeros(Np,size(parameter_new,2));
    parfor k=1:Np    
        % read the parameter_new(k,:) to certain variable
        current_param=parameter_new(k,:);
        parameter_median2=Generatep_fold(parameter_new,k,bound(1,:),bound(2,:));
        alp_old=betas_new*target(current_param);
        alp_media=betas_new*target(parameter_median2);      
 
        % M-H accept
        ratio=min([1,exp(alp_media-alp_old)]);
        u=rand;
        if u<ratio
            result_param=parameter_median2;        % accept
        else
            result_param=current_param; % keep invairable
        end   
        % transfer the result_param to parameter_new(k,:)
        parameter_new_updated(k,:)=result_param;
    end    
    parameter_new=parameter_new_updated;

    % 3.5 save parameter results
    parameter_old=parameter_new;
    parameter_iteration(:,:,stage)=parameter_old;
    disp(['mission has completed :', num2str(stage/S)]);

end

% Parallel settings
delete(gcp('nocreate'));% close parallel pool
end


function [betas]= sequencesGen(S,m,x)
% S         the number of generation 
% Expoential {beta_s} sequence (E. Jeremiah et al., 2012)
s=0:S;
s=s';

% using the relationship: phas=m*s^x to determine x 
betas=m*s.^x;

end

function parameter_new=Generatep_fold(parameter_old,k,LB,UB) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Propose: Generate a candidate particle  for MCMC crossover procedure
%           using the differential evolution Markov Chain algorithm (ter Braak, 2006)
%  Specifically: If the proposal falls outside the range, apply the "fold"
%           operator to fold the proposal back into the valid range
%  References: Jasper A.Vrugt 2016 
[N,d]=size(parameter_old);
gamma=2.38/sqrt(2*d); 
b=1e-4; 

% select two particle that is different current particle k
B=setdiff(1:1:N,k);
R=B(randperm(N-1,2));
parameter_new=parameter_old(k,:)+gamma.*(parameter_old(R(1),:)-parameter_old(R(2),:))+b*randn;

% If out of range, perform "fold" operation
for i=1:d
    % Fold dimensions that are greater than the upper bound
    while parameter_new(i)>UB(i)
        parameter_new(i)=LB(i)+(parameter_new(i)-UB(i));
    end
    % Fold dimensions that are less than the lower bound
    while parameter_new(i)<LB(i)
        parameter_new(i)=UB(i)-(LB(i)-parameter_new(i));
    end
end
end

function outIndex =ResampSys(w, N) 
% Draw a total of N samples with probabilities proportional to the weight vector w, using Systematic Resampling algorithm.
% w            : normalized weight vector (sum to one)
% N (optional) : total number of samples; default to length(w)
% outIndex     : each element is an index into w, or, the "parent" of
%                the sample. Therefore if {X, w} is the original 
%                particles-weights pair, then {X(outIndex), 1/N}
%                will be the resampled pair.   
%
% Author: Lingji Chen
% Date: December 30, 2005
 
eps = 1e-12;            % small but not too small 
len = length(w);
F = cumsum(w);
if abs(F(end) - 1) > eps
  error('the weight vector should be normalized.');
end
 
switch nargin
 case 1
  N = len;
 case 2
 otherwise
  error('wrong number of arguments');
end
 
s = rand / N; 
inc = 1 / N;
 
outIndex = zeros(1, N);
j = 1;
for i = 1:N
  while F(j) < s
    j = j + 1;
  end
  outIndex(i) = j;
  s = s + inc;
end
end

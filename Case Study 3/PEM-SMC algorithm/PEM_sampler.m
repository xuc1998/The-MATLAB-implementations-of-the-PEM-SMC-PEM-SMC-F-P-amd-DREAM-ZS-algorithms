function [paramter_iteration]=PEM_sampler(Np,S,bound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PEM-SMC algorithm used to sample from the target function
% Date: Nov 27, 2021
% Authors: Zhu Gaofeng
% At Lanzhou University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pc=0.6;                             % crossover propabilty
dem=size(bound,2);                  % the dimesion of the target function  
paramter_iteration=nan(Np,dem,S);   % matric used to save samplers duirng the iteration
observed_LE=importdata('./Observations/LE_verify.txt');   % the filepath of the ynthetic observation of LE
observed_NEE=importdata('./Observations/NEE_verify.txt'); % the filepath of the ynthetic observation of NEE
observed_SWC=importdata('./Observations/RSM_verify.txt'); % the filepath of the ynthetic observation of NEE
observed_data=[observed_LE,observed_NEE,observed_SWC];

%% Initialize initial population
% Generate initial population from the sampling space 
for i=1:Np
    for d=1:dem
        parameter_old(i,d)=bound(1,d)+(bound(2,d)-bound(1,d))*rand;
    end 
    w(i,1)=1/Np;
end
 
% Prior distribution; Assuming that the prior distribution of parameters is Uniform 
L_P0=0;
LB=bound(1,:);
UB=bound(2,:); 
for i=1:length(LB)
    L_P0=L_P0-log(UB(i)-LB(i));%?
end
 
% Identify the sequence {beta_s}
[betas]=sequencesGen(S);
betas(S+1:S+5)=1;        % the last five values are forced to be 1


%% SMC Iteration
paramter_iteration(:,:,1)=parameter_old;
for stage=2:S
    % Reweight:                 
    % iteration for Np particle
    betas_new=betas(stage);
    betas_old=betas(stage-1);    
    cal_weight=zeros(Np,1);
    
    for k=1:Np
        target1=target(parameter_old(k,:),observed_data,k);
        alp=betas_new*target1;       % log  of target distribution
        prix=(1-betas_new)*L_P0;     % log  of prior distribution    
        
        alp_old=betas_old*target1;
        prix_old=(1-betas_old)*L_P0;           
        cal_weight(k,1)=(alp+prix)-(alp_old+prix_old);  
    end
    w=w.*exp(cal_weight);
    w=w./(sum(w));
    
   %============================================== Resampling
   Neff(stage,1)=1./sum(power(w,2));
   if Neff(stage,1)<=0.8*Np
       disp('OK');
       % resampling
       ind=ResampSys(w, Np);
       parameter_new=parameter_old(ind,:);    
       w=repmat((1/Np),Np,1); 
       % random walk   
       for i=1:Np
            parameter_median=randw(parameter_new(i,:),LB,UB);
            alp_old=betas_new*target(parameter_new(i,:),observed_data,i);
            alp_media=betas_new*target(parameter_median,observed_data,i);    
             % M-H accept
            ratio=min([1,exp(alp_media-alp_old)]);
            u=rand;
            if u<ratio
                parameter_new(i,:)=parameter_median;  
            else
                parameter_new(i,:)=parameter_new(i,:); 
            end        
       end 
       % crossover operator
       [parameter_median,index]=crossover(parameter_new,pc); 
       for k=1:Np/2        
            % Crossover step
            old_ind=index(2*k-1:2*k); 
            alp_old=betas_new*(target(parameter_new(old_ind(1),:),observed_data,k)+target_normal(parameter_new(old_ind(2),:),observed_data,k));               
            alp_media=betas_new*(target(parameter_median(2*k-1,:),observed_data,k)+target_normal(parameter_median(2*k,:),observed_data,k)); 
            % M-H accept
            ratio=min([1,exp(alp_media-alp_old)]);
            u=rand;
            if u<ratio            
                parameter_new(old_ind(1),:)=parameter_median(2*k-1,:);  
                parameter_new(old_ind(2),:)=parameter_median(2*k,:); 
            else
                parameter_new(old_ind(1),:)=parameter_new(old_ind(1),:); 
                parameter_new(old_ind(2),:)=parameter_new(old_ind(2),:); 
            end
       end
       % Mutation Operater 
       parameter_new_updated=zeros(Np,size(parameter_new,2));
        for k=1:Np    
            % read the parameter_new(k,:) to certain variable
            current_param=parameter_new(k,:);
            parameter_median2=Generatep(parameter_new,k,bound(1,:),bound(2,:));
            alp_old=betas_new*target(current_param,observed_data,k);
            alp_media=betas_new*target(parameter_median2,observed_data,k);      
     
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
   else
       parameter_new=parameter_old;
   end
    disp(['Mission has been complished :' num2str(stage*100/S) ' %']);
    %============================================ Save Parameters
    parameter_old=parameter_new;
    paramter_iteration(:,:,stage)=parameter_old;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [betas]= sequencesGen(S)
% S         the number of generation 
% Method I
% Expoential {beta_s} sequence (E. Jeremiah et al., 2012)
s=0:S;
s=s';
m=7e-16; 

% using the relationship: phas=m*s^x to determine x 
x=log(1/m)/log(S);  
betas=m*s.^x;

% Method II
% s=0:S;
% s=s';
% betas=s./S;
% % Method III
% s=0:S;
% s=s';
% betas=log10(s./S*9+1);
end

function [new_gen,indexx]=crossover(old_popu,pc) 
% Crossover step
[~,indexx]=sort(rand(size(old_popu,1),1));
match_gene=old_popu(indexx,:);
pair=size(match_gene,1)/2;
b=1e-6;
dem=size(match_gene,2);
bits=size(match_gene,2);
capirs=rand(pair,1)<pc;
cpoints=randi(bits,[pair,1]);
cpoints=cpoints.*capirs;
for i=1:pair    
    new_gen([2*i-1,2*i],:)=[match_gene([2*i-1,2*i],1:cpoints(i)),match_gene([2*i,2*i-1],cpoints(i)+1:end)];
end
end 



function parameter_new=Generatep(parameter_old,k,LB,UB) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Propose: generate a candidate partilce for MCMC mutation procedure
%   using the different evolution Markov Chain algorithm (ter Braak, 2006)
%
%%%%%%%%%%%%%%%%%%% Method I Ter Braak
[N,d]=size(parameter_old);
gamma=2.38/sqrt(2*d); 
b=1e-4; 
while 1
    % select two particle that is different current particle k
    B=setdiff(1:1:N,k);
    R=B(randperm(N-1,2));

    parameter_new=parameter_old(k,:)+gamma.*(parameter_old(R(1),:)-parameter_old(R(2),:))+b*randn;
    lbjuest=parameter_new>LB;
    ubjuest=parameter_new<UB;
    if sum(lbjuest)==d && sum(ubjuest)==d 
        break        
    end    
end
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


function parameter_new=randw(parameter_old,LB,UB)
[~,dem]=size(parameter_old);
% Method 1
% b=10^(-2)*ones(1,dem);
% e=-b+2*b.*rand(1,dem);
% method II
% cov=2.38/sqrt(2*dem)*eye(dem);

cov=10^(-4)*eye(dem);
while 1
    % method 1
%     parameter_new=parameter_old+e;
    % method II
    parameter_new=mvnrnd(parameter_old,cov);

    lbjuest=parameter_new>LB;
    ubjuest=parameter_new<UB;
    if sum(lbjuest)==dem && sum(ubjuest)==dem
        break
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

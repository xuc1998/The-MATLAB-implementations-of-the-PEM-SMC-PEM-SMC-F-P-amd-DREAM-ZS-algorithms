function [log_L] = target(x,plugin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: calculating the normal likelihood value of a given parameter
% x for observed data
% x: parameter combination
% plugin: the observations of three targets: LE, NEE, RSM
% log_L is the log transfer of the target distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% observed_data 
observed_data(:,1)=plugin.observed_LE;
observed_data(:,2)=plugin.observed_NEE;
observed_data(:,3)=plugin.observed_RSM;
T=size(observed_data,1);

% file path of the CoLM model execution programs
old_path='/CoLM model execution program';

% use parpool to generate new_path
new_path=[old_path,'-',num2str(randi(1e10))];

% Copy directory 
copyfile(old_path,new_path);

% write the x into certain path txt file
fid=fopen([new_path '/CoLM/input_step.txt'],'w');
for i=1:size(x,2)
    fprintf(fid,'%24.16e\n',x(i));
end
fclose(fid);

% run the CoLM model in the Linux environment
system([new_path '/run'])

% read the model output and caculate likelihood value of target LE
model_data_LE=importdata([new_path '/CoLM/output_LE.txt']);
sum_error_LE=sum((model_data_LE-observed_data(:,1)).^2);
sigma_LE=sum_error_LE/T;
log_L_LE=(-T/2*(log(2*pi*sigma_LE))-sum_error_LE/(2*sigma_LE));

% read the model output and caculate likelihood value of target NEE
model_data_NEE=importdata([new_path '/CoLM/output_NEE.txt']);
sum_error_NEE=sum((model_data_NEE-observed_data(:,2)).^2);
sigma_NEE=sum_error_NEE/T;
log_L_NEE=(-T/2*(log(2*pi*sigma_NEE))-sum_error_NEE/(2*sigma_NEE));

% read the model output and caculate likelihood value of target RSM
model_data_RSM=importdata([new_path '/CoLM/output_RSM.txt']);
sum_error_RSM=sum((model_data_RSM-observed_data(:,2)).^2);
sigma_RSM=sum_error_RSM/T;
log_L_RSM=(-T/2*(log(2*pi*sigma_RSM))-sum_error_RSM/(2*sigma_RSM));

% calculate the  joint likelihood log_L
log_L=log_L_LE+log_L_NEE+log_L_RSM;


% delete new directory
rmdir(new_path,'s');

end
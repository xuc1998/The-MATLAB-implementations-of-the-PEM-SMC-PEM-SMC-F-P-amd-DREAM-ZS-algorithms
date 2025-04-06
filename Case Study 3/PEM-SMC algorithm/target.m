function log_L=target(x,observed_data,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: calculating the normal likelihood value of a given parameter
% x for observed data
% x is the given parameter combinations
% log_L is the log transfer of the target distribution
% k is the index used for creating different directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data size
T=size(observed_data,1);

%  file path of the CoLM model execution programs
old_path='/CoLM model execution program/';
new_path=[old_path,num2str(k)];


% Copy directory 
copyfile(old_path,new_path);

% wirte the x into certain path txt file
fid=fopen([new_path '/coLM/input_step.txt'],'w');
for i=1:size(x,2)
    fprintf(fid,'%24.16e\n',x(i));
end
fclose(fid);

% run the CoLM model in the Linux environment
system([new_path '/run']);


% caculate likelihood value of target LE
model_data_LE=importdata([new_path '/CoLM/output_LE.txt']);
sum_error_LE=sum((model_data_LE-observed_data(:,1)).^2);
sigmav_LE=sum_error_LE/T;
log_L_LE=(-T/2*(log(2*pi*sigmav_LE))-sum_error_LE/(2*sigmav_LE));

% caculate likelihood value of target NEE
model_data_NEE=importdata([new_path '/coLM/output_NEE.txt']);
sum_error_NEE=sum((model_data_NEE-observed_data(:,2)).^2);
sigmav_NEE=sum_error_NEE/T;
log_L_NEE=(-T/2*(log(2*pi*sigmav_NEE))-sum_error_NEE/(2*sigmav_NEE));

% caculate likelihood value of target RSM
model_data_RSM=importdata([new_path '/CoLM/output_RSM.txt']);
sum_error_RSM=sum((model_data_RSM-observed_data(:,3)).^2);
sigmav_RSM=sum_error_RSM/T;
log_L_RSM=(-T/2*(log(2*pi*sigmav_RSM))-sum_error_RSM/(2*sigmav_RSM));

% calculate the  joint likelihood log_L
log_L=log_L_LE+log_L_NEE+log_L_RSM;

%delete new directory
rmdir(new_path,'s');


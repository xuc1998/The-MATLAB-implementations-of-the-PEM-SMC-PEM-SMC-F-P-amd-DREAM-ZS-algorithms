function L=target(x)
%% Description
%  This fucntion is used to calcualte the likelihood  for a given x
%  x is the sample
%  L is the log transfer of the target distribution

%% two-dimensional probability distribution with 20 modes
Nm=20;                        % number of mode   
sigma=ones(1,20)*0.1^2;       % variance of normal distribution
mu=[2.18 5.76;8.67 9.59;4.24 8.48;8.41 1.68;3.93 7.82;3.25 3.47;1.70 0.50;
    4.59 5.60;6.91 5.81;6.87 5.40;5.41 2.65;2.70 7.88;4.98 3.70;1.14 2.39;
    8.33 9.50;4.93 1.50;1.83 0.09;2.26 0.31;5.54 6.86;1.69 8.11];
w=ones(1,Nm)*0.05;            % weights
f=0;
for i=1:Nm
    f=f+w(i)/(2*pi*sigma(i))*exp(-(x-mu(i,:))*(x-mu(i,:))'/(2*sigma(i)));
end
L=log(f);

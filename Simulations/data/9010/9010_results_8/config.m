%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% config.m: Configuration file of the INEX tool
% author: Kerstin Lenk
% date: 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory for the resulting spike trains
string_directory = directory;

% Length of the spike train in [msec]
lengthST = 300000;

% time interval in [sec] including refractory time
t = 0.005;         

% Number of excitatory neurons
ex0 = 200; 
ex1 = 0;

% Number of inhibitory neurons
in0 = 50;
in1 = 0;

% Connectivity (between 0 and 100, decimal number)
%conn = 20;

% Lower and upper boundary of random number for ... 
% ... basic activity
c_k = [0.02;0.0;0.02]; 

maxSynStr = 0.7;
% ... synaptic strength(all positive numbers)
synStr_oex0 = [maxSynStr;0.2;maxSynStr]; % excitatory, to other neurons, group 0
synStr_oex1 = [0.0;0.0;0.0]; % excitatory, to other neurons, group 1
synStr_sex0 = [0.0;0.0;0.0]; % excitatory, autapses, group 0
synStr_sex1 = [0.0;0.0;0.0]; % excitatory, autapses, group 1
synStr_oin0 = [maxSynStr;0.2;maxSynStr]; % inhibitory, to other neurons, group 0
synStr_oin1 = [0.0;0.0;0.0]; % inhibitory, to other neurons, group 1
synStr_sin0 = [0.0;0.0;0.0]; % inhibitory, autapses, group 0
synStr_sin1 = [0.0;0.0;0.0]; % inhibitory, autapses, group 1

sensitivityMultiplier = 0.1;







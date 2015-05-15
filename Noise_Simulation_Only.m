% Noise Simulation Only is the top-level simulation code that implements 
% and executes the noise-simulating algorithms and prepares the PSF's for
% final analysis and PSF subtraction.

% Author: Josh Samuels
% March 2015



clear all;
close all;
%cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'
% Limit Comp Threads for parrallel processing
maxNumCompThreads(72);
load('TSM_Config.mat');

total_pics = M_num_exposures + M_num_reference;
modifier = floor(M_num_exposures/M_num_reference);
m = 1;
num_exp = 0;
num_ref = 0;
while m <= total_pics 
    for n = 1:modifier
        if (num_exp < M_num_exposures)
            Noise_Simulation(m,true);
            m = m+1;
            num_exp = num_exp + 1;
        end
    end
    if (num_ref < M_num_reference)
        Noise_Simulation(m,false);
        m = m+1;
        num_ref = num_ref + 1;
    end
end

% load('TSM_Config.mat');
% for m = 1:M_num_exposures
%     Noise_Simulation(m,true);
% end
% 
% for m = 1:M_num_reference
%     Noise_Simulation((m+M_num_exposures),false);
% end


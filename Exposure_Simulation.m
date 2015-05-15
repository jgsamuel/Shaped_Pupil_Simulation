% Exposure Simulation is the top-level simulation code that implements and
% executes the physical modeling of CCD exposures using a space-based 
% telescope with a coronagraph system.

% Exposure Simulation can be used to add together the PSF's to 
% and generate a simulation of the entire target system, generating sample 
% data as if it was collected by a simulated camera.

% Author: Josh Samuels
% March 2015
% Object positions should be defined in terms of diffraction widths
% (units of lambda/D)



clear all;
close all;
%cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'
load('TSM_Config.mat');

total_pics = M_num_exposures + M_num_reference;
modifier = floor(M_num_exposures/M_num_reference);
m = 1;
num_exp = 0;
num_ref = 0;
while m <= total_pics 
    for n = 1:modifier
        if (num_exp < M_num_exposures)
            Intensity_Modeling_Only('TSM_Config.mat',m);
            m = m+1;
            num_exp = num_exp + 1;
        end
    end
    if (num_ref < M_num_reference)
        Reference_Star_Only('TSM_Config.mat', m);
        m = m+1;
        num_ref = num_ref + 1;
    end
end

% % Simple pattern of exposures, not interspersed
% Generate Science Exposures
% for m = 1:M_num_exposures
%    Intensity_Modeling_Only('TSM_Config.mat',m);
% end
% 
% % Generate Reference Star Exposures
% for m = 1:M_num_reference
%     Reference_Star_Only('TSM_Config.mat', (m+M_num_exposures));
% end

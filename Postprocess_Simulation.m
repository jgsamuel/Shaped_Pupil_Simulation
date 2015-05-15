% Postprocess Simulation is the top-level simulation code that implements 
% and executes the reading and compilation of data exposures, and then 
% performs analysis and PSF subtraction on the resulting data.
% 

% Author: Josh Samuels
% March 2015



clear all;
close all;
%cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'
% Limit Comp Threads for parrallel processing
maxNumCompThreads(72);
load('TSM_Config.mat');
load('Psn_Sim_Out_1.mat');

% Cutout populated section of CCD. 
[row, col] = find(CCD_Expectation);
edge = sqrt(length(col));
%r_size = size(CCD_Expectation,1);
%c_size = size(CCD_Expectation,2);

Science_Target_Frames = zeros(edge,edge,M_num_exposures);
Reference_Star_Frames = zeros(edge,edge,M_num_reference);
Average_Max = 0;
Total_Frames = M_num_exposures + M_num_reference;

modifier = floor(M_num_exposures/M_num_reference);
m = 1;
num_exp = 0;
num_ref = 0;
while m <= Total_Frames 
    for n = 1:modifier
        if (num_exp < M_num_exposures)
            inname = (['Psn_Sim_Out_' num2str(m) '.mat']);
            load(inname);
            
            Science_Target_Frames(:,:,(num_exp+1)) = Random_Data(min(row):max(row),min(col):max(col)); % Replace with Random_Data
            Average_Max = max(max(Random_Data))+Average_Max;
            clearvars 'CCD_Expectation' 'Random_Data'
            
            m = m+1;
            num_exp = num_exp + 1;
        end
    end
    if (num_ref < M_num_reference)
        inname = (['Ref_Sim_Out_' num2str(m) '.mat']);
        load(inname);
        
        Reference_Star_Frames(:,:,(num_ref+1)) = Random_Data(min(row):max(row),min(col):max(col)); % Replace with Random_Data
        Average_Max = max(max(Random_Data))+Average_Max;
        clearvars 'CCD_Expectation' 'Random_Data'
 
        m = m+1;
        num_ref = num_ref + 1;
    end
end

% for m = 1:M_num_exposures
%     inname = (['Psn_Sim_Out_' num2str(m) '.mat']);
%     load(inname);
%     %Science_Target_Frames(:,:,m) = CCD_Expectation(min(row):max(row),min(col):max(col)); % Replace with Random_Data
%     Science_Target_Frames(:,:,m) = Random_Data(min(row):max(row),min(col):max(col)); % Replace with Random_Data
%     Average_Max = max(max(Random_Data))+Average_Max;
%     clearvars 'CCD_Expectation' 'Random_Data'
% end
% 
% for m = 1:M_num_reference
%     inname = (['Ref_Sim_Out_' num2str(m + M_num_exposures) '.mat']);
%     load(inname);
%     %Reference_Star_Frames(:,:,m) = CCD_Expectation(min(row):max(row),min(col):max(col)); % Replace with Random_Data
%     Reference_Star_Frames(:,:,m) = Random_Data(min(row):max(row),min(col):max(col)); % Replace with Random_Data
%     Average_Max = max(max(Random_Data))+Average_Max;
%     clearvars 'CCD_Expectation' 'Random_Data'
% end

%Cleaned_Science_Target = zeros(edge,edge,M_num_exposures);
[Cleaned_Science_Target, Z, sigma] = PCA_PSF_Subtraction('TSM_Config.mat', Reference_Star_Frames, Science_Target_Frames, 10);
Dumb_PSF_Subtraction_Frames = Cleaned_Science_Target;
Reference_Frame = squeeze(sum(Reference_Star_Frames,3))./(size(Reference_Star_Frames,3));
num_frames = size(Science_Target_Frames,3);
for f = 1:num_frames
    Dumb_PSF_Subtraction_Frames(:,:,f) = Dumb_PSF_Subtraction_Frames(:,:,f)-Reference_Frame;
end
Dumb_PSF_Image = squeeze(sum(Dumb_PSF_Subtraction_Frames,3));

% Add together PSF-subtracted science data
Science_Image = squeeze(sum(Cleaned_Science_Target,3));

Stacked_Image = squeeze(sum(Science_Target_Frames,3));
Stacked_Image = Stacked_Image./(size(Science_Target_Frames,3));

Average_Max = Average_Max/Total_Frames;
save('Reduced_Data.mat', 'Dumb_PSF_Image', 'Stacked_Image', 'Science_Image', 'Cleaned_Science_Target', 'Z', 'sigma', 'Average_Max');




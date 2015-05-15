
clear all;
close all;
%cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'
% Limit Comp Threads for parrallel processing
maxNumCompThreads(72);
load('TSM_Config.mat');
load('Psn_Sim_Out_1.mat');
clearvars 'pupil_plane_aberration';

% Cutout populated section of CCD. 
[row, col] = find(CCD_Expectation);
edge = sqrt(length(col));
%r_size = size(CCD_Expectation,1);
%c_size = size(CCD_Expectation,2);
%M_num_exposures = 70;
%M_num_reference = 20;

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
            
            Science_Target_Frames(:,:,(num_exp+1)) = CCD_Expectation(min(row):max(row),min(col):max(col)); % Replace with Random_Data
            Average_Max = max(max(Random_Data))+Average_Max;
            clearvars 'CCD_Expectation' 'Random_Data'
            
            m = m+1;
            num_exp = num_exp + 1;
        end
    end
    if (num_ref < M_num_reference)
        inname = (['Ref_Sim_Out_' num2str(m) '.mat']);
        load(inname);
        
        Reference_Star_Frames(:,:,(num_ref+1)) = CCD_Expectation(min(row):max(row),min(col):max(col)); % Replace with Random_Data
        Average_Max = max(max(Random_Data))+Average_Max;
        clearvars 'CCD_Expectation' 'Random_Data'
 
        m = m+1;
        num_ref = num_ref + 1;
    end
end

Average_Max = Average_Max/Total_Frames;
nFrames = size(Science_Target_Frames,3);

videoName = 'Dynamic_Noise_Sim_w_Psn';
%filepath = 'C:\Users\jsamuels\Documents\MATLAB\Thesis\Animations\';


imageMat = cell(1,nFrames);
for p = 1:nFrames
    
end


save('Vid_Temp_Out.mat','-v7.3');
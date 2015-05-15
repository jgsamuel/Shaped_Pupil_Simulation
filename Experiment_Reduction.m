% Experiment Reduction is a working script for analyzing K-variation in the PCA-subtracted data
clear all;
load('Psn_Sim_Out_Speckles_17K.mat');
load('PCA_ExperimentLambda17_2_2_200.mat');
load('Reduced_Data_Speckles_17K.mat');
% Annulus Detection

% Planet Location:
r_Planet = 99; % Row location of target object
c_Planet = 110; % Column location of target object
num_rows = size(Science_Image,1);
num_cols = size(Science_Image,2);
center_row = num_rows./2;
center_col = num_cols./2;

r_p = r_Planet - center_row;
r_c = c_Planet - center_col;
annulus_radius = r_p.*r_p + r_c.*r_c;
max_a_r = ceil(sqrt(annulus_radius)+1.5);
min_a_r = floor(sqrt(annulus_radius)-1.5);

count = 1;
k_vec = 2:2:100;
pca_vec = k_vec;
for k = 2:2:100
    k
    Science_Image = Science_Images(:,:,count);
    
    Planet_Intensity = Science_Image(r_Planet, c_Planet);

    annulus = [];

    for r = 1:num_rows
        for c = 1:num_cols
            % Determine if pixel is in annulus:
            r_centered = r - center_row;
            c_centered = c - center_col;
            rad = sqrt(r_centered*r_centered + c_centered*c_centered);
            if (rad >= min_a_r && rad <= max_a_r)
                % Pixel in radius
                annulus = [annulus Science_Image(r, c)];
                %disp(['Row: ' num2str(r) '  Col: ' num2str(c)]);
                %Science_Image(r,c) = 10000000; % Draws circle with annulus
            end
        end
    end
    mean_annulus_intensity = mean(annulus);
    stddev_annulus_intensity = std(annulus);
    Observation = (Planet_Intensity-mean_annulus_intensity)/stddev_annulus_intensity
    pca_vec(count) = Observation;
    count = count+1;
end
% figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, Science_Image, [3500 4500]);
% colorbar; axis square;
% %xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('Millimeters');
% title('Annulus Imposed on Binned Science Image');
% set(gca, 'Ydir', 'normal');
% colormap gray;
plot(k_vec, pca_vec)
title('Tuning the Number of Principal Components')
xlabel('K: Number of Principal Components')
ylabel('Signal to Noise Ratio of Earth Observation')
axis square;

annulus_stacked = [];

Planet_Intensity_stacked = Stacked_Image(r_Planet, c_Planet).*size(Cleaned_Science_Target,3);
for r = 1:num_rows
    for c = 1:num_cols
        % Determine if pixel is in annulus:
        r_centered = r - center_row;
        c_centered = c - center_col;
        rad = sqrt(r_centered*r_centered + c_centered*c_centered);
        if (rad >= min_a_r && rad <= max_a_r)
            % Pixel in radius
            annulus_stacked = [annulus_stacked (Stacked_Image(r, c).*size(Cleaned_Science_Target,3))];
            %disp(['Row: ' num2str(r) '  Col: ' num2str(c)]);
            %Science_Image(r,c) = 10000000; % Draws circle with annulus
        end
    end
end

mean_annulus_intensity_stacked = mean(annulus_stacked);
stddev_annulus_intensity_stacked = std(annulus_stacked);
Stacked_Observation = (Planet_Intensity_stacked-mean_annulus_intensity_stacked)/stddev_annulus_intensity_stacked


% fitswrite(Principal_Components,'PC.fits');
% fitswrite(Science_Image,'Science_Image.fits');
% fitswrite(Cleaned_Science_Target, 'PSF_Subtracted_Data.fits');
%  Code wraps a loop around data reduction, finds optimal annulus for data
%  processing meta-analysis:
% width_vec = 0.1:0.1:2.5;
% PCA_Vec = 0.1:0.1:2.5;
% Stacked_Vec = 0.1:0.1:2.5;
% Subtracted_Vec = 0.1:0.1:2.5;
% COUNTER = 1;
% for WIDTH = 0.1:0.1:2.5
%     WIDTH 



% Data Reduction is a working script for analyzing PCA-subtracted data
a = Science_Image./(Average_Max*size(Cleaned_Science_Target,3));
b = Stacked_Image/(Average_Max);
c = (Dumb_PSF_Image-mean(mean(Dumb_PSF_Image)))./(Average_Max*size(Cleaned_Science_Target,3));
a(a<= 0) = eps;
b(b<= 0) = eps;
c(c<= 0) = eps;

close all

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, log10(a), [-10 -9]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('Stacked Science Images with PCA-Enhanced Subtraction, Log-Scale Contrast');
set(gca, 'Ydir', 'normal');
%colormap gray;

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, log10(b), [-10 -9]); % Also try -11 to -9.5
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('Stacked Science Images, Log-Scale Contrast');
set(gca, 'Ydir', 'normal');
%colormap gray;

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, log10(c), [-10 -9]); % Also try -11 to -9.5
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('Stacked Science Images with PSF Subtraction, Log-Scale Contrast');
set(gca, 'Ydir', 'normal');
%colormap gray;

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, (Stacked_Image.*size(Cleaned_Science_Target,3)), [0 16000]); % Also try -11 to -9.5
%figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, log10(c), [-10 -9]); % Also try -11 to -9.5
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('Stacked Science Images, Linear Scale Flux');
set(gca, 'Ydir', 'normal');
colormap gray;

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, Dumb_PSF_Image, [-8000 8000]); % Also try -11 to -9.5
%figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, log10(c), [-10 -9]); % Also try -11 to -9.5
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('Stacked Science Images with Simple PSF Subtraction, Linear Scale Flux');
set(gca, 'Ydir', 'normal');
colormap gray;

figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, Science_Image, [-8000 8000]);% Good stretch: -5K to 4.5K
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('Stacked Science Images with PCA-Enhanced PSF Subtraction, Linear Scale Flux ');
set(gca, 'Ydir', 'normal');
colormap gray;

Principal_Components = zeros(size(Science_Image,1),size(Science_Image,2),size(Z,1));
height = size(Science_Image,1);
width = size(Science_Image,2);
for im = 1:size(Z,1)
    for r = 1:size(Science_Image,1)
        Principal_Components(r,:,im) = Z(im,(((r-1)*width)+1):(r*width));
    end
end


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
max_a_r = ceil(sqrt(annulus_radius)+1.5);%+WIDTH);
min_a_r = floor(sqrt(annulus_radius)-1.5);%-WIDTH);


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
            %Science_Image(r,c) = 10000000; % DrawsW circle with annulus
        end
    end
end

% figure('visible', 'on'); imagesc(x_cam_phys, y_cam_phys, Science_Image, [3500 4500]);
% axis square;
% %xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('Millimeters');
% ylabel('Millimeters');
% title('Annulus Imposed on Binned Science Image');
% set(gca, 'Ydir', 'normal');
% colormap gray;


mean_annulus_intensity = mean(annulus);
stddev_annulus_intensity = std(annulus);
Observation = (Planet_Intensity-mean_annulus_intensity)/stddev_annulus_intensity

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


annulus_subtracted = [];
Planet_Intensity_subtracted = Dumb_PSF_Image(r_Planet, c_Planet);
for r = 1:num_rows
    for c = 1:num_cols
        % Determine if pixel is in annulus:
        r_centered = r - center_row;
        c_centered = c - center_col;
        rad = sqrt(r_centered*r_centered + c_centered*c_centered);
        if (rad >= min_a_r && rad <= max_a_r)
            % Pixel in radius
            annulus_subtracted = [annulus_subtracted Dumb_PSF_Image(r, c)];
            %disp(['Row: ' num2str(r) '  Col: ' num2str(c)]);
            %Science_Image(r,c) = 10000000; % Draws circle with annulus
        end
    end
end

mean_annulus_intensity_subtracted = mean(annulus_subtracted);
stddev_annulus_intensity_subtracted = std(annulus_subtracted);
Subtracted_Observation = (Planet_Intensity_subtracted-mean_annulus_intensity_subtracted)/stddev_annulus_intensity_subtracted


% fitswrite(Principal_Components,'PC.fits');
% fitswrite(Science_Image,'Science_Image.fits');
% fitswrite(Cleaned_Science_Target, 'PSF_Subtracted_Data.fits');


% Code to test optimal annulus for data post-analysis:
% PCA_Vec(COUNTER) = Observation;
% Stacked_Vec(COUNTER) = Stacked_Observation;
% Subtracted_Vec(COUNTER) = Subtracted_Observation;
% COUNTER = COUNTER + 1;
% end % end for loop for trying different annulus widths
% close all
% hold on
% plot(width_vec, PCA_Vec)
% plot(width_vec, Stacked_Vec)
% plot(width_vec, Subtracted_Vec)
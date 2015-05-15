% % ViewOutput extracts and displays the output of the coronagraph and
% error simulation code. It reads Psn_Sim_Out_x files (x = exposure
% number), and plots the output at various stages of the simulation.

% Input: Exposure number. Enter 0 to view 'Psn_Sim_Out.mat'

% Author: Josh Samuels
% March 2015
function [success] = ViewOutput(EXPOSURE_NUMBER)
close all;
cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'
if (EXPOSURE_NUMBER == 0)
    infile = 'Psn_Sim_Out.mat';
else
    infile = (['Psn_Sim_Out_' num2str(EXPOSURE_NUMBER) '.mat'])
end
if (exist(infile) == 0)
    infile = (['Ref_Sim_Out_' num2str(EXPOSURE_NUMBER) '.mat'])
    if (exist(infile) == 0)
        success = false;
        return;
    end
end
success = true;
load(infile);
% Plot Solar System PSF
fig = figure('visible', 'off');
plot(X_Cut_Vec, X_Cut_Intensity);
ylim([-11 0]);
title('Solar System , band-averaged');
xlabel('lambda/D Radians');
savefig('Center_Row_Contrast');

fig = figure('visible', 'off'); imagesc(XVEC_CUT, XVEC_CUT, log10(IMG_PLANE_INTENSITY), [-11 -9.5]);
colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('lambda/D Radians');
title('Solar System (Angle on the sky), band-averaged');
set(gca, 'Ydir', 'normal');
savefig('Solar_System_Bandavg_Angle')

% Plot physical units version
fig = figure('visible','off'); imagesc(xvec_phys, xvec_phys, log10(IMG_PLANE_INTENSITY), [-11 -9.5]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('Solar System Image Plane (physical units), band-averaged');
set(gca, 'Ydir', 'normal');
savefig('Solar_System_Bandavg_Phys')

% % Plot physical units version
% fig = figure('visible','off'); imagesc(xvec_phys, xvec_phys, log10(IMG_PLANE_SAMP), [-11 -9.5]);
% colorbar; axis square;
% %xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('Millimeters');
% title('Solar System w/ Speckle Noise (physical units), band-averaged');
% set(gca, 'Ydir', 'normal');
% savefig('Solar_System_Rician_Sample')

% % Plot physical units version of smoothed Rician noise 
% fig = figure('visible','off'); imagesc(xvec_phys, xvec_phys, log10(IMG_PLANE_SMOOTHED), [-11 -9.5]);
% colorbar; axis square;
% %xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('Millimeters');
% title('Solar System w/ Smoothed Speckle Noise (physical units), band-averaged');
% set(gca, 'Ydir', 'normal');
% savefig('Solar_System_Rician_Smoothed')


figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(CCD_Expectation), [0 1.5]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('CCD Expectation(physical units)');
set(gca, 'Ydir', 'normal');
savefig('CCD_Expectation');

% Plot Random Data matrix

% Remove negative realizations
Random_Data(Random_Data <= 0) = eps;
Random_Data = Random_Data./(max(max(Random_Data)));
fig = figure('visible','off'); imagesc(x_cam_phys, y_cam_phys, log10(Random_Data), [-10.5 -9]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
title('Solar System Data Sample w/ Noise Factors, band-averaged');
set(gca, 'Ydir', 'normal');
savefig('Random_Data')

% Random_Samp = poissrnd(CCD_Expectation.*1e13);
% Random_Samp = Random_Samp./(max(max(Random_Samp))); %normalize
% figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(Random_Samp), [-10.5 -9]);
% colorbar; axis square;
% %xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('Millimeters');
% title('CCD Sample(physical units)');
% set(gca, 'Ydir', 'normal');
% savefig('Random_Sample');


close all;

% View mask
figure; imagesc(gray2d_SP); colorbar; colormap gray; axis square; title('binned-down gray SP');

% Show figures
openfig('Center_Row_Contrast', 'new', 'visible');
openfig('Solar_System_Bandavg_Angle', 'new', 'visible');
openfig('Solar_System_Bandavg_Phys', 'new', 'visible');
%openfig('Solar_System_Rician_Sample', 'new', 'visible');
openfig('Random_Data', 'new', 'visible');
openfig('CCD_Expectation', 'new', 'visible');
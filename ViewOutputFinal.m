% % ViewOutput extracts and displays the output of the coronagraph and
% error simulation code. It reads Psn_Sim_Out_x files (x = exposure
% number), and plots the output at various stages of the simulation.

% Input: Exposure number. Enter 0 to view 'Psn_Sim_Out.mat'

% Author: Josh Samuels
% March 2015
close all;
% Phase Screen Plot:
% imagesc(1:2000, 1:2000, PHASE_SCREEN, [-5e-4 5e-4]);
% colorbar; axis square;
% xlabel('Pupil Plane Pixels');
% ylabel('Pupil Plane Pixels');
% title('Phase Screen Intensity');
% set(gca, 'Ydir', 'normal');

% Plot Solar System PSF
fig = figure('visible', 'off');
plot(X_Cut_Vec, X_Cut_Intensity);
ylim([-11 0]);
title('Solar System , band-averaged');
xlabel('lambda/D Radians');
title('Log-Scale Image Plane Contrast (Normalized to Peak)')
ylabel('Relative Intensity')
savefig('Center_Row_Contrast');

fig = figure('visible', 'off'); imagesc(XVEC_CUT, XVEC_CUT, log10(IMG_PLANE_INTENSITY), [-10.5 -8.5]);
colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('lambda/D Radians');
ylabel('lambda/D Radians');
title('Solar System (Angle on the sky), band-averaged');
set(gca, 'Ydir', 'normal');
savefig('Solar_System_Bandavg_Angle')

% Plot physical units version
fig = figure('visible','off'); imagesc(xvec_phys, xvec_phys, log10(IMG_PLANE_INTENSITY), [-10.5 -8.5]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
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

figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(CCD_Expectation/(max(max(CCD_Expectation)))), [-10.5 -8.5]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('CCD Expectation (units of contrast)');
set(gca, 'Ydir', 'normal');
savefig('CCD_Expectation');

figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, CCD_Expectation, [0 150]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('CCD Expectation (photon count)');
set(gca, 'Ydir', 'normal');
savefig('CCD_Expectation_Flux');

% Plot Random Data matrix

% Remove negative realizations
Random_Data(Random_Data <= 0) = eps;
Random_Data = Random_Data./(max(max(Random_Data)));
fig = figure('visible','off'); imagesc(x_cam_phys, y_cam_phys, log10(Random_Data), [-10.5 -8.5]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('Solar System Data Sample w/ Noise Factors, band-averaged');
set(gca, 'Ydir', 'normal');
savefig('Random_Data')

Super_Sample = poissrnd(CCD_Expectation.*10);
Super_Sample = Super_Sample./(max(max(Super_Sample))); %normalize
figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(Super_Sample), [-11 -9.5]);%[-10.5 -9]);
colorbar; axis square;
%xlim([-OPA OPA]); ylim([-OPA OPA]);
xlabel('Millimeters');
ylabel('Millimeters');
title('CCD Sample(Oversaturated)');
set(gca, 'Ydir', 'normal');
savefig('Super Exposed Image');


close all;

% View mask
figure; imagesc(gray2d_SP); colorbar; colormap gray; axis square; title('binned-down gray SP');

% Show figures
openfig('Center_Row_Contrast', 'new', 'visible');
openfig('Solar_System_Bandavg_Angle', 'new', 'visible');
openfig('Solar_System_Bandavg_Phys', 'new', 'visible');
openfig('CCD_Expectation', 'new', 'visible');
openfig('CCD_Expectation_Flux', 'new', 'visible');
openfig('Super Exposed Image', 'new', 'visible');
openfig('Random_Data', 'new', 'visible');

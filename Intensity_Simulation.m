% % Intensity_Simulation simulates the PSF of a single celestial object using
% the input mask. Its output is intended to represent the expectation of 
% photon count at any position in the image plane.
%
% This code is mostly adapted from code developed in November 2014 by 
% Neil Zimmerman in order to compute polychromatic PSF's for shaped pupil 
% coronagraphs (eval_gray2d_crmask_polychromPSF)
% 
% % Author: Josh Samuels
% March 2015

function [Target_Object_Normalized_PSF, CORE_PSF_THROUGHPUT, CORE_PSF, TOTAL_PSF_INTENSITY] = Intensity_Simulation(MASK_FILENAME, PHASE_SCREEN, OUTER_WORKING_ANGLE, DIFFRACTION_WIDTHS_X, DIFFRACTION_WIDTHS_Y, SAMPLE_FAC, FRAC_BW)

% Configure for text output
outfileID = fopen('output.txt','w');

%% Load mask
Target_Object_Normalized_PSF = NaN;
CORE_PSF_THROUGHPUT = NaN;
TOTAL_PSF_INTENSITY = NaN;
IWA = 6.;
%OWA = 60;   % Outer working angle
OWA = OUTER_WORKING_ANGLE; % Test to see larger space
OPA = OWA + 10; % outer propagation angle

%ovsampfac = 1;
ovsampfac = SAMPLE_FAC;
M_pup = 1000;
L = 2*ovsampfac*(2*M_pup);

fracBW = FRAC_BW;
% Nlambda = 31;
Nlambda = 7;
gamma = linspace(1-fracBW/2, 1+fracBW/2, Nlambda);

gray2d_SP_fitsname = MASK_FILENAME; %sprintf('TPF_CRMASK_6WA60_C10_localopt_gray_D%04d.fits', 2*M_pup);
gray2d_SP = fitsread(gray2d_SP_fitsname);
% figure; imagesc(gray2d_SP); colorbar; colormap gray; axis square; title('binned-down gray SP');

%% central wavelength only
% on-axis
%
E_in_onax = ones(2*M_pup).*exp(1i*PHASE_SCREEN);
E_clean_onax = ones(2*M_pup);
E_SP_onax = gray2d_SP.*E_in_onax;
E_clean_SP_onax = gray2d_SP.*E_clean_onax;
E_SP_pad_onax = zeros(L);
E_SP_pad_onax(L/2-(M_pup-1):L/2+M_pup, L/2-(M_pup-1):L/2+M_pup) = E_SP_onax;
% Record pupil plane graphs
%E_SP_pad_onax_clean = zeros(L);
%E_SP_pad_onax_clean(L/2-(M_pup-1):L/2+M_pup, L/2-(M_pup-1):L/2+M_pup) = E_clean_SP_onax;
%I_pup_onax = abs(E_SP_pad_onax).^2;
%I_pup_onax_clean = abs(E_SP_pad_onax_clean).^2;
%save('Pupil_Plane_Pics.mat','I_pup_onax','I_pup_onax_clean');
%imagesc(I_pup_onax)

E_foc_onax = 1/L*fftshift( fft2( fftshift(E_SP_pad_onax) ) );
E_foc_onax_peak = sum(sum(E_SP_onax))*(1/(2*ovsampfac))*(0.5/M_pup);

% off-axis
%
[X_pup,Y_pup] = meshgrid(-(M_pup-1/2):1:(M_pup-1/2), -(M_pup-1/2):1:(M_pup-1/2));
% tilt = 10.; % diffraction widths
%tilt = DIFFRACTION_WIDTHS_X;
tilt_x = DIFFRACTION_WIDTHS_X; % diffraction widths (lambda/D)
tilt_y = DIFFRACTION_WIDTHS_Y;

% offax_phase_ramp = 2*pi*X_pup/(2*M_pup)*tilt;
% x_offax = L/2+1 + 2*ovsampfac*tilt;
% y_offax = L/2+1;
D_x = X_pup/(2*M_pup);
D_y = Y_pup/(2*M_pup);

offax_phase_ramp_x = 2*pi*D_x*tilt_x;
offax_phase_ramp_tot = offax_phase_ramp_x + 2*pi*D_y*tilt_y;

E_in_offax = exp(1i*(offax_phase_ramp_tot + PHASE_SCREEN));
E_SP_offax = gray2d_SP.*E_in_offax;
E_SP_pad_offax = zeros(L);
E_SP_pad_offax(L/2-(M_pup-1):L/2+M_pup, L/2-(M_pup-1):L/2+M_pup) = E_SP_offax;

E_foc_offax = 1/L*fftshift( fft2( fftshift(E_SP_pad_offax) ) );

%% Make plots of central wavelength PSF
%
I_foc_onax = abs(E_foc_onax).^2;
% I_foc_masked_onax = abs(E_foc_masked_onax).^2;

I_foc_offax = abs(E_foc_offax).^2;
% I_foc_masked_offax = abs(E_foc_masked_offax).^2;

img_cutout_beg = L/2+1 - 2*ovsampfac*OPA;
img_cutout_end = L/2+1 + 2*ovsampfac*OPA;
center_row = 2*ovsampfac*OPA + 1;
I_foc_onax_cutout = I_foc_onax(img_cutout_beg:img_cutout_end, img_cutout_beg:img_cutout_end);
I_foc_offax_cutout = I_foc_offax(img_cutout_beg:img_cutout_end, img_cutout_beg:img_cutout_end);

% figure; imagesc(log10(I_foc_onax_cutout/max(max(I_foc_onax_cutout))), [-9 0]); colorbar; axis square;
xvec_cut = (-2*ovsampfac*OPA:2*ovsampfac*OPA)/(2*ovsampfac);
[XX_cut, YY_cut] = meshgrid(xvec_cut);

%Added work to find core of PSF. Note that if PSF is not rotationally
%symmetric, the boundary finding code below will probably not work.
%Instead, delete everything up to IS_IN_CORE and just use logical indexing 
% on the entire intensity matrix.
Center_Row_Axis = XX_cut(center_row, center_row:end);
FWHM_Boundary = 0; % Mean with half max of PSF, determine "core" in terms of lambda/D in main wavelength
FWHM_Index = 1;
Intensity_Radius = I_foc_onax_cutout(center_row, center_row:end)/(E_foc_onax_peak^2);
for r = 1:length(Center_Row_Axis)
    if (Intensity_Radius(r) > 0.5)
        FWHM_Boundary = Center_Row_Axis(r);
        FWHM_Index = r;
    end
end
%disp(['Full Width Half-Max at ' num2str(MWHM_Boundary) ' lambda/D'])
outtext = ['Full Width Half-Max at ' num2str(FWHM_Boundary) ' lambda/D'];
fprintf(outfileID, '%s\n',outtext);
%CORE_PSF = I_foc_onax_cutout(center_row, (center_row + MWHM_Index));
core_cutout_beg = L/2+1 - 2*ovsampfac*FWHM_Boundary;
core_cutout_end = L/2+1 + 2*ovsampfac*FWHM_Boundary;
CORE_PSF = I_foc_onax(core_cutout_beg:core_cutout_end, core_cutout_beg:core_cutout_end);
IS_IN_CORE = ((CORE_PSF./(max(max(I_foc_onax)))) < 0.5);  
CORE_PSF(IS_IN_CORE) = 0;
TOTAL_PSF_INTENSITY = sum(sum(I_foc_onax));
CORE_PSF_THROUGHPUT = sum(sum(CORE_PSF))/TOTAL_PSF_INTENSITY;
%disp(['Core PSF Throughput: ' num2str(CORE_PSF_THROUGHPUT)])
outtext = ['Core PSF Throughput: ' num2str(CORE_PSF_THROUGHPUT)];
fprintf(outfileID, '%s\n',outtext);

% figure; plot(XX_cut(center_row, center_row:end), log10(I_foc_onax_cutout(center_row, center_row:end)/E_foc_onax_peak^2)); ylim([-11 0]);
% title('Contrast at first focal plane, central wavelength');
% xlabel('lambda/D');

%off axis PSF
% figure; imagesc(xvec_cut, xvec_cut, log10(I_foc_offax_cutout), [-3 0]); colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('lambda/D');
% title('Off-axis PSF, central wavelength');

%% compute polychromatic PSF
%
clearvars 'E_SP_pad_onax' 'E_SP_pad_offax' 'E_foc_offax' 'E_foc_onax'

I_foc_polychrom_onax_cutout = zeros(length(gamma), 2*ovsampfac*OPA*2 + 1, 2*ovsampfac*OPA*2 + 1);
I_foc_polychrom_offax_cutout = zeros(length(gamma), 2*ovsampfac*OPA*2 + 1, 2*ovsampfac*OPA*2 + 1);
I_foc_polychrom_onax_peak = zeros(length(gamma),1);
for gii = 1:length(gamma),
    lambda_ratio = gamma(gii);
    sL = round(ceil(L*lambda_ratio/2)*2);

    E_in_onax = ones(2*M_pup).*exp(1i*PHASE_SCREEN/lambda_ratio);
    E_SP_onax = gray2d_SP.*E_in_onax;
    E_SP_pad_onax = zeros(sL);
    E_SP_pad_onax(sL/2-(M_pup-1):sL/2+M_pup, sL/2-(M_pup-1):sL/2+M_pup) = E_SP_onax;
    E_foc_onax_peak = sum(sum(gray2d_SP))/sL;
    I_foc_polychrom_onax_peak(gii) = E_foc_onax_peak^2;

    E_SP_pad_offax = zeros(sL);
    %offax_phase_ramp = 2*pi*X_pup/(2*M_pup) * tilt/lambda_ratio;
    offax_phase_ramp_x = 2*pi*D_x * tilt_x;
    offax_phase_ramp_tot = offax_phase_ramp_x + 2*pi*D_y * tilt_y;
    total_phase_shift = offax_phase_ramp_tot + PHASE_SCREEN;
    E_in_offax = exp(1i*(total_phase_shift/lambda_ratio));
    E_SP_offax = gray2d_SP.*E_in_offax;
    E_SP_pad_offax(sL/2-(M_pup-1):sL/2+M_pup, sL/2-(M_pup-1):sL/2+M_pup) = E_SP_offax;
    
    E_foc_onax = 1/sL*fftshift( fft2( fftshift(E_SP_pad_onax) ) );
    I_foc_onax = abs(E_foc_onax).^2;
    E_foc_offax = 1/sL*fftshift( fft2( fftshift(E_SP_pad_offax) ) );
    I_foc_offax = abs(E_foc_offax).^2;
    
    img_cutout_beg = sL/2+1 - 2*ovsampfac*OPA;
    img_cutout_end = sL/2+1 + 2*ovsampfac*OPA;
    center_row = 2*ovsampfac*OPA + 1;
    I_foc_polychrom_onax_cutout(gii,:,:) = I_foc_onax(img_cutout_beg:img_cutout_end, img_cutout_beg:img_cutout_end);
    I_foc_polychrom_offax_cutout(gii,:,:) = I_foc_offax(img_cutout_beg:img_cutout_end, img_cutout_beg:img_cutout_end);
    clearvars 'E_SP_pad_onax' 'E_SP_pad_offax' 'E_foc_offax' 'E_foc_onax'
end

xvec_cut = (-2*ovsampfac*OPA:2*ovsampfac*OPA)/(2*ovsampfac);
[XX_cut, YY_cut] = meshgrid(xvec_cut);

%% Make plots of band-averaged PSF
% Plot radial contrast
I_foc_bandavg_onax_cutout = squeeze(mean(I_foc_polychrom_onax_cutout, 1));
I_foc_bandavg_offax_cutout = squeeze(mean(I_foc_polychrom_offax_cutout, 1));
I_foc_bandavg_onax_peak = squeeze(mean(I_foc_polychrom_onax_peak));
%I_foc_bandavg_offax_peak = squeeze(mean(I_foc_polychrom_offax_peak)); Well
%approximated with onax

I_foc_bandavg_norm_onax = I_foc_bandavg_onax_cutout / I_foc_bandavg_onax_peak;
I_foc_bandavg_norm_offax = I_foc_bandavg_offax_cutout / I_foc_bandavg_onax_peak;

Target_Object_Normalized_PSF = I_foc_bandavg_norm_offax;


% 
% figure;
% plot(XX_cut(center_row, center_row:end), log10(I_foc_bandavg_norm_onax(center_row, center_row:end)));
% ylim([-11 0]);
% title('Contrast in image plane, band-averaged');
% xlabel('lambda/D');
% 
% % Plot on axis, 2-D PSF
% figure; imagesc(xvec_cut, xvec_cut, log10(I_foc_bandavg_norm_onax), [-11 0]);
% colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('lambda/D');
% title('Normalized, ideal, band-averaged star PSF');
% 
% % Plot Sun + Earth at maximum elongation
% I_sun_plus_earth = I_foc_bandavg_norm_onax + I_foc_bandavg_norm_offax/1e10;
% 
% plot(XX_cut(center_row, center_row:end), log10(I_sun_plus_earth(center_row, center_row:end)));
% ylim([-11 0]);
% title('Sun + Earth, band-averaged');
% xlabel('lambda/D');
% 
% figure; imagesc(xvec_cut, xvec_cut, log10(I_sun_plus_earth), [-11 -9.5]);
% colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
% xlabel('lambda/D');
% title('Sun + Earth, band-averaged');



fclose(outfileID);
end


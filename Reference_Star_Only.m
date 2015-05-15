% % Reference_Star_Only is a modification of Intensity_Modeling_Only, which
% computes the expected PSF's of a list of celestial bodies through the input mask.

% Instead, however, this function will create images without any planetary
% bodies, simulating the collection of "reference star" data where planets
% are ignored or smoothed away. While in reality there are multiple ways of
% creating such photographs by smoothing away 

% WARNING: RAM usage scales extremely quickly in terms of "ovsampfac".
% Running at a decent resolution (e.g. 3-5) requires 10-40 GB of RAM. Use
% Exposure_Simulation to implement this code on a 

% Note: Object positions should be defined in terms of diffraction 

% Author: Josh Samuels
% March 2015


function [outname] = Reference_Star_Only(configFile,EXPOSURE_NUMBER)
    
load(configFile);
    
    phase_screen = pupil_plane_aberration(:,:,EXPOSURE_NUMBER);
    % Alternatively, generate simple phase screen:
    %phase_screen = ABB2D((2*M_pup), 140, 200, 3, (2*pi/1000), true); %(2*pi/16000), true);
    %phase_screen = Generate_Phase_Screens(1, (2*M_pup), OPA, 200, (2*pi/2000), (2*pi/16000));
    % Or, no phase screen at all:
    %phase_screen = zeros(2*M_pup);
    
    [Image_Plane_Intensity, T_PSF, CORE_PSF, TOTAL_PSF_INTENSITY] = Intensity_Simulation(gray2d_SP_fitsname, phase_screen, OWA, positionsX(1),positionsY(1), ovsampfac, fracBW);

    % Naive approach: Need to sum over the entire area of the PSF, cutout
    % is not negligible. 
    %Center_Point_Fraction = Image_Plane_Intensity(center_row, center_row)/(sum(sum(Image_Plane_Intensity)));
    Center_Point_Fraction = max(max(CORE_PSF))/TOTAL_PSF_INTENSITY;
    
    % Compute Intensity from Starlight
    % All PSF's are normalized to the max point of the star PSF. Star flux is given
    % by space simulation in photons/sec. De-normalize by multiplying by
    % flux at that point. 
    F_star = photon_count*Center_Point_Fraction*exposure_time;
    
%     No object PSF's in reference image
%     for o = 2:objects
%         I_mat = Intensity_Simulation(gray2d_SP_fitsname, phase_screen, OWA, positionsX(o), positionsY(o), ovsampfac, fracBW);
%         Image_Plane_Intensity = Image_Plane_Intensity + I_mat * magnitudes(o);
%     end

    clearvars 'I_mat'

    X_Cut_Vec = XX_cut(center_row, center_row:end);
    X_Cut_Intensity = log10(Image_Plane_Intensity(center_row, center_row:end));
    XVEC_CUT = xvec_cut;
    IMG_PLANE_INTENSITY = Image_Plane_Intensity;
    % % Plot Solar System PSF
    % fig = figure('visible', 'off');
    % plot(XX_cut(center_row, center_row:end), log10(Image_Plane_Intensity(center_row, center_row:end)));
    % ylim([-11 0]);
    % title('Solar System , band-averaged');
    % xlabel('lambda/D Radians');
    % savefig('Center_Row_Contrast');
    % 
    % fig = figure('visible', 'off'); imagesc(xvec_cut, xvec_cut, log10(Image_Plane_Intensity), [-11 -9.5]);
    % colorbar; axis square; xlim([-OPA OPA]); ylim([-OPA OPA]);
    % xlabel('lambda/D Radians');
    % title('Solar System (Angle on the sky), band-averaged');
    % set(gca, 'Ydir', 'normal');
    % savefig('Solar_System_Bandavg_Angle')

    % Compute Intensity from Starlight: Old method
    % All PSF's are normalized to the flux of the star. Star flux is given
    % by space simulation in photons/sec, 
    % F_star = photon_count*exposure_time;
    
    % Add Intensity from Exozodi Profile
    % Assume uniform distribution of zodi within the working area of the star
    % F_ez = irradiance density of flux, derived in Space Simulation function
    % Method integrates over all zodi, using Parseval's theorem to sum up over the
    % entire PSF. Result can be multiplied by area of each pixel
    F_ez = F_zodi; % See Kasdin's mean photon count, equivalent contrast, and integration time due to exo-zodi as a function of telescope diameter (unpublished as of Spring 2015)

    A_psuedoarea = (sum(sum(gray2d_SP)))/(size(gray2d_SP,1)*size(gray2d_SP,1)); %Normalized psuedoarea
    A_psuedoarea = A_psuedoarea * pi * ((Telescope_Diameter/2)^2);  % Multiply psuedoarea by area of actual telescope
    %Incorrect with normalized PSF, based on Kasdin eq. 97:
    %Intensity_Exozodi = F_ez * max(max(Image_Plane_Intensity)) * (lambda*lambda) * (1/(A_psuedoarea));
    Intensity_Exozodi = F_ez*(1/(focal_length*focal_length)) * (A_psuedoarea);
    Intensity_Exozodi = Intensity_Exozodi*exposure_time;
    % Save workspace, as the following part has low space requirements and high
    % runtime requirements that can be optimized using parralellization
    clearvars 'pupil_plane_aberration'
    outname = (['REF_TEMP_OUT_' num2str(EXPOSURE_NUMBER) '.mat']);
    save(outname);
end
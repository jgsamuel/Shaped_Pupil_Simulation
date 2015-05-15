% % This file runs the second half of the simulation, the parts of
% which are runtime-heavy and can be optimized using parrallelization

% TODO:  All
%    
% Object positions should be defined in terms of diffraction widths
% (units of lambda/D)


% Note: If is_Science = false, the function assumes the exposure is a
% reference image and looks for the file accordingly.
function [outname] = Noise_Simulation(EXPOSURE_NUMBER, is_Science)
    
    %clear all;
    if (is_Science)
        inname = (['SIM_TEMP_OUT_' num2str(EXPOSURE_NUMBER) '.mat']);
        load(inname); % Load the completed space-heavy component of MATLAB simulation
    else
        inname = (['REF_TEMP_OUT_' num2str(EXPOSURE_NUMBER) '.mat']);
        load(inname); % Load the completed space-heavy component of MATLAB simulation
    end
    % 20 threads maximum on Adroit Cluster. Restriction is disliked by MATLAB but is imperative for
    % good behavior while on cluster.
    %maxNumCompThreads(20);
    %pixelScale = 20*um;
    
    % Flux calculation from Space Simulation (moved to Intensity_Simulation
    % code)
    %F_star = photon_count*exposure_time;

    % Add speckle noise from a Rician Speckle Noise distribution, as discussed
    % in Soummer's 2007 paper. Warning: Extremely runtime-intensive.
    % This produces an even plane of noise, however, as there is no theoretical
    % definition of spatial correlation between pixels.
    % Image_Plane_Intensity = Rician_Speckle_Noise_Simulation(Image_Plane_Intensity, OWA, ovsampfac, fracBW);
    % IMG_PLANE_SAMP = Image_Plane_Intensity;
    % IMG_PLANE_SMOOTHED = PSF_Convolution(IMG_PLANE_SAMP, CORE_PSF);
    IMG_PLANE_SAMP = Image_Plane_Intensity;
    IMG_PLANE_SMOOTHED = Image_Plane_Intensity;

    % Define expected photon counts at each pixel of the CCD

    % Vanderbei Fun Suggestion:
    % Camera problems? Flats, darks.  CCD cameras have a finite well depth.
    % When they fill, they spill over into the neighboring pixel. Saturated
    % pixel will actually make a spike into the overflowing pixels. Pixel
    % overflow should have a simple pattern, fall down (in direction that the
    % system reads. 

    % Reformulate PSF into physical units.
    % lambda/D radians (angle on the sky) is equivalent to (lambda * f)/ D in
    % physical distance across the image plane

    xvec_phys = xvec_cut*(lambda).*(focal_length).*(1/Telescope_Diameter);
    xvec_phys = xvec_phys * 1000; % convert to mm
    % fig = figure('visible','off'); imagesc(xvec_phys, xvec_phys, log10(Image_Plane_Intensity), [-11 -9.5]);
    % colorbar; axis square;
    % %xlim([-OPA OPA]); ylim([-OPA OPA]);
    % xlabel('Millimeters');
    % title('Solar System Image Plane (physical units), band-averaged');
    % set(gca, 'Ydir', 'normal');
    % savefig('Solar_System_Bandavg_Phys')


    %CCD_Expectation = imresize(Image_Plane_Intensity,(1/ovsampfac));
    % Map pixels to Image_Plane_Intensity positions
    CCD_Expectation = zeros(ccdHeightPix,ccdWidthPix);
    Tracker = zeros(ccdHeightPix,ccdWidthPix);
    x_cam = ((-ccdHeightPix/2):1:(ccdHeightPix/2)); % x axis in pixels (rows)
    y_cam = ((-ccdWidthPix/2):1:(ccdWidthPix/2)); % y axis in pixels (cols)
    pixel_mm = pixelScale * 1000; % pixel in millimeters
    x_cam_phys = x_cam.*pixel_mm;
    y_cam_phys = y_cam.*pixel_mm;

    
    % Convert Image_Plane_Intensity to a Flux Matrix
    Image_Plane_Flux = Image_Plane_Intensity.*F_star;
    Pixels_Total = 0;
    Pixels_N = 0;
    % Integrate over intensity profile for pixel (x,y):
    % which is in between x_cam_phys(x) and x_cam_phys(x+1) and between 
    % y_cam_phys(y) and y_cam_phys(y+1)
    current_physical_position_x = NaN;
    current_physical_position_y = NaN;
    current_physical_position_xPlusOne = x_cam_phys(1);
    current_physical_position_yPlusOne = y_cam_phys(1);
    image_index_x_start = NaN;
    image_index_x_end = 1; 
    image_index_y_start = NaN;
    image_index_y_end = 1;
    for x = 1:ccdHeightPix
        % Update physical location into a start and ending index
        current_physical_position_x = current_physical_position_xPlusOne;
        current_physical_position_xPlusOne = x_cam_phys(x+1);
        image_index_x_start = image_index_x_end; % Maybe add one, so that no overlap? But when sampling is low probably a bad idea
        if (current_physical_position_x > xvec_phys(length(xvec_phys)))
            continue; % Ignore pixels outside of the PSF
        end

        if (sum(xvec_phys < current_physical_position_xPlusOne) == 0)
            continue;% Ignore pixels outside of the PSF
        end
        image_index_x_end = sum(xvec_phys < current_physical_position_xPlusOne);
        image_index_y_end = 1;
        for y = 1: ccdWidthPix
            % Update Y axis physical 
            current_physical_position_y = current_physical_position_yPlusOne;
            current_physical_position_yPlusOne = y_cam_phys(y+1);
            image_index_y_start = image_index_y_end; % Maybe add one, so that no overlap? But when sampling is low probably a bad idea
            if (current_physical_position_y > xvec_phys(length(xvec_phys)))
                continue; %Ignore pixels outside of the PSF
            end
            if (sum(xvec_phys < current_physical_position_yPlusOne) == 0)
                continue; %Ignore pixels outside of the PSF
            end
            image_index_y_end = sum(xvec_phys < current_physical_position_yPlusOne);

            % Integrate over intensity region
            Pixels_Sum = 0;
            Pixels_Total = Pixels_Total+1;
            Pixel_Components = 0;
            for a = image_index_x_start:image_index_x_end
                for b = image_index_y_start:image_index_y_end
                    Pixels_N = Pixels_N + 1;
                    Pixel_Components = Pixel_Components + 1;
                    Pixels_Sum = Pixels_Sum + Image_Plane_Flux(a,b);
                end
            end
            
            %CCD_Expectation(x,y) = CCD_Expectation(x,y) + (Pixels_Sum/Pixel_Components) + Intensity_Exozodi*Pixel_Area; 
            CCD_Expectation(x,y) = (Pixels_Sum/Pixel_Components);
            Tracker(x,y) = Pixel_Components;
        end
    end

    [row, col] = find(CCD_Expectation);
    edge = sqrt(length(col));
    Pixel_Area = pixelScale*pixelScale;
    
    % Figure out average number of simulated pixels in each physical pixel
    Sim_Pixel_Per_CCD_Pixel = Pixels_N/Pixels_Total;
    
    CCD_Expectation = CCD_Expectation(min(row):max(row),min(col):max(col));  % Cut out pixels outside of the simulation area
    x_cam_phys = x_cam_phys(min(row):max(row));
    y_cam_phys = y_cam_phys(min(col):max(col));
    
    CCD_Expectation = CCD_Expectation.*Sim_Pixel_Per_CCD_Pixel; % Scale up flux for binning effect
    CCD_Expectation = CCD_Expectation+Intensity_Exozodi*Pixel_Area; % Add Per-Pixel Flux for Exozodi
    
    % figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(CCD_Expectation), [-11 -9.5]);
    % colorbar; axis square;
    % %xlim([-OPA OPA]); ylim([-OPA OPA]);
    % xlabel('Millimeters');
    % title('CCD Expectation(physical units)');
    % set(gca, 'Ydir', 'normal');
    % savefig('CCD_Expectation');

    % Take sample photograph
    %images = 53;
    rows = size(CCD_Expectation,1);
    col = size(CCD_Expectation,2);
    DARK_FRAME_FLAT_FLUX = 100;
    Mean_Dark_Frame = DARK_FRAME_FLAT_FLUX.*(ones(rows,col));
    READ_NOISE_FLAT_STD_DEV = 1; % Electron multiplication CCD have extremely low read noise
    Read_Noise_Var = READ_NOISE_FLAT_STD_DEV.*(ones(rows,col));
    Read_Noise_Mean = zeros(rows,col);
    
    % Dark frame addition: Add gain factor? X is psn process of photon
    % arrival, but actually reads Y = gamma*X.   Can determine this because
    % variance will be gamma^2*X 
    
    % Generate Simulated Image
    dark_sample = poissrnd(Mean_Dark_Frame)-Mean_Dark_Frame;
    read_noise = normrnd(Read_Noise_Mean, Read_Noise_Var);
    Random_Data = poissrnd(CCD_Expectation) + dark_sample + read_noise;
    % Code to write "Image":
    %filename = [files_location 'SampleImage_' num2str(im) '.fits'];
    %fitswrite(Random_Data, filename);

    % Random_Samp = poissrnd(CCD_Expectation.*1e12);
    % Random_Samp = Random_Samp./(max(max(Random_Samp))); %normalize
    % figure('visible', 'off'); imagesc(x_cam_phys, y_cam_phys, log10(Random_Samp), [-11 -9.5]);
    % colorbar; axis square;
    % %xlim([-OPA OPA]); ylim([-OPA OPA]);
    % xlabel('Millimeters');
    % title('CCD Sample(physical units)');
    % set(gca, 'Ydir', 'normal');
    % savefig('Random_Sample');

    if (is_Science)
        outname = (['Psn_Sim_Out_' num2str(EXPOSURE_NUMBER) '.mat']);
        save(outname,'gray2d_SP', 'X_Cut_Vec', 'X_Cut_Intensity', 'Random_Data', 'CORE_PSF', 'XVEC_CUT', 'IMG_PLANE_INTENSITY','IMG_PLANE_SAMP', 'IMG_PLANE_SMOOTHED', 'OPA', 'xvec_phys','x_cam_phys', 'y_cam_phys', 'CCD_Expectation');
    else
        outname = (['Ref_Sim_Out_' num2str(EXPOSURE_NUMBER) '.mat']);
        save(outname,'gray2d_SP', 'X_Cut_Vec', 'X_Cut_Intensity', 'Random_Data', 'CORE_PSF', 'XVEC_CUT', 'IMG_PLANE_INTENSITY','IMG_PLANE_SAMP', 'IMG_PLANE_SMOOTHED', 'OPA', 'xvec_phys','x_cam_phys', 'y_cam_phys', 'CCD_Expectation');
    end
end
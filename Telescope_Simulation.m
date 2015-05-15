% Telescope Simulation is the top-level simulation code that declares the
% simulation parameters, draws random noise, and generates photographs of 
% the celestial model as noise patterns change through time.

% Exposure Simulation can then be used to add together the PSF's to 
% and generate a simulation of the entire target system, generating sample 
% data as if it was collected by a simulated camera.

% Author: Josh Samuels
% March 2015


clear all;
close all;
cd 'C:\Users\jsamuels\Documents\MATLAB\Thesis'

% Limit Comp Threads for parrallel processing
maxNumCompThreads(20);

% Camera parameters  
um = 1e-6;
%pixelScale = 24*um; % Forrestal science camera (2013)
%pixelScale = 8*um; % Canon DSLR t3i
pixelScale = 5*um;
%pixelScale = 1.5*um; %iPhone pixels, sample size
ccdWidthPix = 1024; 
ccdHeightPix = 1024;  


% Experiment Imaging Parameters
M_num_exposures = 30;
M_num_reference = 10;
M_total_pics = M_num_exposures + M_num_reference;
exposure_time = 600; % seconds


% Model and image parameters
IWA = 6.; OWA = 60; 
%OWA = 200; % Test larger area
OPA = OWA + 10; % outer propagation angle

ovsampfac = 1;
M_pup = 1000;
L = 2*ovsampfac*(2*M_pup);

lambda = 550e-9; %meters
focal_length = 100; %meters
Telescope_Diameter = 8; %meters
fracBW = 0.1;   % +/- 5 %
%fracBW = 0.01; 
% Nlambda = 31;
Nlambda = 7;


% Begin formulating input data
img_cutout_beg = L/2+1 - 2*ovsampfac*OPA;
img_cutout_end = L/2+1 + 2*ovsampfac*OPA;
center_row = 2*ovsampfac*OPA + 1;



gamma = linspace(1-fracBW/2, 1+fracBW/2, Nlambda);

gray2d_SP_fitsname = sprintf('TPF_CRMASK_6WA60_C10_localopt_gray_D%04d.fits', 2*M_pup);
gray2d_SP = fitsread(gray2d_SP_fitsname);

%LOCATION_TO_SAVE_PICTURES
files_location = 'C:\Users\jsamuels\Documents\MATLAB\Thesis\2-5\';



xvec_cut = (-2*ovsampfac*OPA:2*ovsampfac*OPA)/(2*ovsampfac);
[XX_cut, YY_cut] = meshgrid(xvec_cut);


[photon_count, flux_ratio, positionsX, positionsY, F_zodi] = Space_Simulation(exposure_time, M_num_exposures, gray2d_SP_fitsname, lambda, focal_length, Telescope_Diameter, fracBW);
% magnitudes = [1 (1/1e10) (1/1e10)];
% positionsX = [0 10 46.33];
% positionsY = [0 10 0];
magnitudes = [1 flux_ratio];
positionsX = [0 positionsX];
positionsY = [0 positionsY];

objects = length(magnitudes);

% Generate Random Time-Dependent Aberrations

num_screens = 12;
%
%phase_screens = zeros((2*M_pup),(2*M_pup),num_screens); 
%phase_screen = ABB2D((2*M_pup), 60, 200, 3, (2*pi/2000), true); %(2*pi/16000), true);
%phase_screens = Generate_Phase_Screens(num_screens, (2*M_pup), OWA, 100,(2*pi/4000), (2*pi/34000)); %Original numbers
phase_screens = Generate_Phase_Screens(num_screens, (2*M_pup), OWA, 100, (2*pi/200), (2*pi/34000)); % Based on polishing law

% calculate total time spent observing a target for time scaling
TOTAL_TIME = M_total_pics * exposure_time; 
functions = cell(num_screens,1);

% Make some errors increase monotonically over time
num_increasing = floor(num_screens/4);
functions{1} = @(x) (x/(TOTAL_TIME))^2;
for f = 2:num_increasing
    modifier = f-3;
    functions{f} = @(x) (x/(TOTAL_TIME/(2^(modifier))))-(2^(modifier-1));
end

num_decreasing = floor(num_screens/4);
functions{num_increasing+1} = @(x) ((-1)*((x/(TOTAL_TIME))^2));
for f = (num_increasing+2):(num_increasing+num_decreasing)
    modifier = f-(num_increasing+3);
    functions{f} = @(x) (2^(modifier-1))-(x/(TOTAL_TIME/(2^(modifier))));
end

% Make the rest of the functions periodic
% sin((2*pi*t)/T) oscillates with period T

% BASE PERIOD (seconds)
period = exposure_time*24.3; % Not even fraction of exposure times

for f = (num_increasing+num_decreasing+1):num_screens
    functions{f} = @(x) sin((2*pi*x)/period);
    period = period*2;
end

% Configure for exposures:
time = 0; % total seconds elapsed
%weights = ones(M_num_exposures,num_screens);
%weights = ones(num_screens,1);
pupil_plane_aberration = zeros((2*M_pup),(2*M_pup),M_num_exposures);
%total_W = []; % track weights through time
for m = 1:M_total_pics
    full_screen = zeros((2*M_pup),(2*M_pup));
    %w = []; %Track weights through time
    for s = 1:num_screens
        weight = arrayfun(functions{s},time);
        %w = [w weight];
        to_add = weight.*phase_screens(:,:,s);
        full_screen = full_screen + to_add;
    end
    %w       % Track weights through time
    %sum(abs(w))
    %total_W = [total_W sum(abs(w))];

    pupil_plane_aberration(:,:,m) = full_screen;
    time = time + exposure_time;
end
clearvars 'full_screen' 'phase_screens' 'to_add' 'functions'
save('TSM_Config.mat', '-v7.3');
%  Run Photo Simulations (commented out to allow for parrallelization of
%  Telescope_Simulation.m)
% for m = 1:M_num_exposures
%     %Intensity_Modeling_Only('TSM_Config.mat',m);
% end

% % Space_Simulation processes astronomical data into a unitless format
% that can be used for coronagraph modeling. 
function [photon_count, flux_ratio, positionsX, positionsY, F_ez] = Space_Simulation(EXPOSURE_TIME, NUM_EXPOSURES, MASK_FILENAME, lambda, focal_length, Telescope_Diameter,FRAC_BW)

% Parameters
gray2d_SP = fitsread(MASK_FILENAME);
%Telescope_Diameter = 8; %meters
%lambda = 550e-9; %meters
fracBW = FRAC_BW;


% Planet Information (Source: HORIZONS http://ssd.jpl.nasa.gov/horizons.cgi
% Mercury:  Geometric albedo is 0.106  Mean Radius is 2440 km
% Venus:  Geometric albedo is 0.65  Mean Radius is 6051.8 km
% Earth: Geometric albedo is 0.367  Mean Radius is 6371.01 km
% Mars: Geometric albedo is 0.150   Mean Radius is 3389.9 km
% Jupiter: Geometric albedo is 0.52   Mean Radius is 69911 km

% Distances to the Sun (also has eccentricity on site), mean distances 
%(elliptical orbits)
% Source: http://hyperphysics.phy-astr.gsu.edu/hbase/Solar/soldata2.html#c5
% Mercury:  0.387 AU
% Venus:  0.723 AU
% Earth: 1 AU  (by definition)
% Mars: 1.524 AU
% Jupiter: 5.203 AU


%r = Seperation between object and Star 
%D = distance to Star
%theta = angle between star and object, as seen from earth

%angular_position = [(pi/2) (pi/2) (pi/2) (pi/2) (pi/2)]; 
%angular_position = [0 0 0 0 0]; 
angular_position = [1.5 1.57 5 5.2 6]; % Arbitrary numbers selected for useful image
%angular_position = [GET EMPHERIS DATA];  % Actual orbits based on data above
%angle within the orbit. Set angle 
%("theta") to be such that theta of zero is the planet being directly
% "behind" the star, when orbit is less than face-on (inclination, or phi, 
% is greater than zero).
% It appears that my definition should be analagous to the "true anomaly"
% http://en.wikipedia.org/wiki/Longitude_of_the_ascending_node

geometric_albedo = [0.106 0.65 0.367 0.150 0.52]; 
%orbital_inclination_deg = [0 0 0 0 0];
orbital_inclination_deg = [66 66 66 66 66];
orbital_inclination = (orbital_inclination_deg./360).*(2*pi); % angle between face on and actual orbital plane
% 0 degrees is face on, 90 degrees is edge-on.  Units in RADIANS.
% Positive inclination means that while "above" the star, from the point
% of view of the observer, it is also "behind" the star, and at lower phase
% angle.  But this doesn't matter too much, as the rotational position of 
% the observer is arbitrary.



r = [0.387 0.723 1 1.524 5.203];    % Actual  physical orbital distance, AU
object_radius = [2440 6051.8 6371.01 3389.9 69911]; % Radius of Planets, Km
object_radius = object_radius.*1000; % Convert to meters
distance_to_object = 7.61;
D = linspace(distance_to_object,distance_to_object, length(r)); % Distance to object, parsecs

% First attempt, true only at periastron: = r.*(cos(orbital_inclination));
% Project r into plane of sky, from point of view of observer
r_x = -sin(angular_position).*r;
r_y = (cos(orbital_inclination)).*(cos(angular_position)).*r;


% 1 AU 
parsecInAU = 206264.806; %Source : Google.com
AUInMeters = 149597870700; % Source: Google.com
% theta ~= tan(theta)  = r/D
theta = r./(D.*parsecInAU);
theta_x = r_x ./ (D.*parsecInAU);
theta_y = r_y ./ (D.*parsecInAU);
% 1 arcsecond = 1/3600 degrees  convert to radians (*(2pi/360)) =
% pi/648000
theta = theta ./ (pi/648000); % Degrees to ArcSeconds
theta_x = theta_x ./ (pi/648000); % Degrees to ArcSeconds
theta_y = theta_y ./ (pi/648000); % Degrees to ArcSeconds
%disp('Milliarcseconds:')
theta = theta * 1000;
theta_x = theta_x * 1000;
theta_y = theta_y * 1000;
lambdaOverD = lambda/Telescope_Diameter;

% Approximate Angular resolution: 
%disp('Milliarcseconds:')
R = 1000* (lambdaOverD ./ (pi/648000));

positions = theta./R;
% Default: Set positions to be in a line, to the right of star
positionsX = theta_x./R;
positionsY = theta_y./R;

% Energy and Flux Calculations
% Source: http://www.neilzim.net/wp-content/uploads/2014/12/AFTA_read_noise.pdf
f_lambda = 3.63*10^-8; % Based off of Johnson-Cousins V-band zero
um = 1e-6;
lambda = 550e-9; %meters
photon_energy = (6.626e-34)*(3e8)*(1/lambda); % hc/lambda 
lambda_micrometers = lambda * (1/um);

absolute_magnitude_star = 4.83; %M_V  Absolute magnitude of Sol.  Source: www.google.com
%distance_to_object = 10; % Already defined above
apparent_mag = absolute_magnitude_star + 5*log10(distance_to_object) - 5;

photon_flux = (f_lambda*lambda_micrometers*fracBW)/photon_energy;
photon_flux = photon_flux*(10^((-1/2.5)*apparent_mag));

quantum_effeciency = 0.8;
eta = 0.97^16; % Estimate of optical effeciency for WFIRST

effective_area = (sum(sum(gray2d_SP)))/(size(gray2d_SP,1)*size(gray2d_SP,1));
effective_area = effective_area * pi * ((Telescope_Diameter/2)^2);
%effective_area = 4.112; % Number used in Neil's sample calculation

% T_psf =Fraction of non-scattered starlight in the core of the PSF
%T_psf = 0.05; % Neil's sample calculation
T_psf = 0.41517;  % Full width half max, calculated at ovsampfac = 6
photon_count = photon_flux*effective_area*eta*quantum_effeciency;
% Using Neil's calculation directly for total photon count over psf:
%photon_count = photon_flux*effective_area*eta*T_psf*quantum_effeciency;
% Using Kasdin's calculation for P(0,0):
%photon_count = photon_flux*effective_area*effective_area*eta*quantum_effeciency*(1/(lambda*lambda))*(1/(focal_length*focal_length));

% Exposure Details
% M_num_exposures = 100;
% exposure_time = 300; % seconds
M_num_exposures = NUM_EXPOSURES;
exposure_time = EXPOSURE_TIME; % seconds
photon_count_tot = photon_count*exposure_time*M_num_exposures;


% Phase Angle= observer planet star orientation, with angle at the planet
% vertex.

% 
% Ellipse transformation stuff: http://web.ncf.ca/aa456/scale/ellipse.html
% http://www.newton.dep.anl.gov/askasci/math99/math99114.htm

phase_angle = acos(cos(angular_position).*sin(orbital_inclination_deg));

phi = sin(phase_angle) + (pi - phase_angle).*cos(phase_angle);
phi = phi./pi;
planet_ratio = object_radius./(r.*AUInMeters);
flux_ratio = geometric_albedo.*(planet_ratio.*planet_ratio);
flux_ratio = flux_ratio.*phi;




% Calculate Flux from Zodi and Exozodi
absolute_magnitude_reference_star = 0; %User zero-mag star for scaling, see Kasdin (6)
distance = 7.61; %parsecs
apparent_mag = absolute_magnitude_reference_star + 5*log10(distance) - 5;
% Identical flux calculation to above
photon_flux = (f_lambda*lambda_micrometers*fracBW)/photon_energy;
photon_zero_point = (f_lambda)/photon_energy; % Without band?
photon_flux = photon_flux*(10^((-1/2.5)*apparent_mag));
photon_zero_point = photon_zero_point*(10^((-1/2.5)*apparent_mag)); % If preserving spectral intensity

Mag_zodi = 23; % Mag/arcsec^2 for exozodi at 1 AU seperation for Earth system 
Mag_exozodi = 22;
exozodi_number = 1; % Number of exozodi in units of local zodi

zodi_flux = 10^((-1/2.5)*Mag_zodi);
zodi_flux = zodi_flux*4.25*10^10; % Convert from arcsec^2 to steradians
exo_flux = 10^((-1/2.5)*Mag_exozodi);
exo_flux = exo_flux*4.25*10^10; % Convert from arcsec^2 to steradians

F_0 = photon_flux; % Assuming band can be fixed before integrating over exozodi on Kasdin p6
F_ez = F_0*(zodi_flux + exozodi_number*exo_flux);
%F_ez = F_ez*eta*T_psf*quantum_effeciency; % Drop effective area from flux
%calc because it will be multiplied by pixel area. PSF throughput T_PSF likely not relevant here.
F_ez = F_ez*eta*quantum_effeciency;
end

% ABB2D generates phase errors that match a Kolmogorov (or any power law) distribution
% Originally written by Alexis Carlotti
% I have modified the power law to fit polishing errors.

% Parameters:
% 1. Width of the output array
% 2. Minimum spatial frequency (cycles per aperture)
% 3. Maximum spatial frequency
% 4. Number of passes through the generator (e.g. 3)
% 5. Scaling factor applied to the screen. In the example below, the factor
%    sets the root-mean-square phase to be 2*pi/20 radians.
%    (Modify to adjust the strength of the phase errors)

% 6. Sets power law to Kolmogorov (true or false)
% Sample call:
%  kolmogorov_screen = ABB2D(2000, 1, 60, 3, 2*pi/20, 1);

function P=ABB2D(N,N1Waves,N2Waves,Nt,RMS,Kolmogorov)

% Perturb2D(N,NWaves,Nt,RMS) returns an N by N matrix whose values have a root mean square of RMS and an offset of 1.
% NWaves is the maximum number of aperture cycles that is considered
% (usually 60)
% Nt is the number of random waves for 1 aperture cycle (usually 3)

%N=500;
x=(-N/2+0.5:N/2-0.5)/N;
[X,Y]=meshgrid(x);

%NWaves=60;
%Nt=3;
%RMS=1;

Phase=zeros(N);
Temp=zeros(N);

for k=N1Waves:1/Nt:N2Waves
    
        alpha=2*pi*rand();
        Temp=Temp+cos(2*pi*((Y+X*tan(alpha))/(sqrt(1+tan(alpha)^2)))*k+2*pi*rand());
    
    if Kolmogorov==1
        % Kolmogorov Power Law
        %Phase=Phase+k^(-5/3)*Temp;
        % Power Law for Optical Polishing Errors
        Phase=Phase+k^(-2)*Temp;
    else
        
        Phase=Phase+Temp;
        
    end
    
end

P=RMS*Phase/sqrt(sum(sum((Phase-mean(mean(Phase))).^2))/N/N);
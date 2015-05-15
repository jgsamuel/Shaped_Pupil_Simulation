%Generate phase screens creates a set of phase screens that can be used to
% simulate a group of phase abberrations that vary with time

% Physical System Parameters
% Actuators: 140 x 140 ( because OPA = 70) (because you need two samples per sinusoidal period to correct for) Nyquist theorem
% 70 ... 200
% phase ripple in frequency is a shifted impulse function in fourier space
% speckle's distance from star is its frequency

function [phase_screens] = Generate_Phase_Screens(num_screens, size, OWA, Outside_Angle, Scaling_Normal, Scaling_Suppressed)

    % Passes through the Kolmogorov generator
    PASSES = 3;
    phase_screens = zeros(size,size,num_screens);

    parfor i = 1:num_screens
        % Add two regimes together--one suppressed by AO, one not. Useful in theory but causes residual speckles
        % in lower spatial frequencies that can be problematic.
        %screen = ABB2D(size, OWA, Outside_Angle, PASSES, Scaling_Normal, true);
        %screen = screen + ABB2D(size, 1, OWA, PASSES, Scaling_Suppressed, true);
        screen = ABB2D(size, 1, OWA, PASSES, Scaling_Suppressed, true);
        phase_screens(:,:,i) = screen;
    end

end
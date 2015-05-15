
clear all;
close all;
load('Vid_Temp_Out_Lambda_3K.mat')
%load('Vid_Temp_Out_Speckles_17K.mat');

videoName = 'Speckles_Lambda3K';
%videoName = 'Speckles_Lambda17K';
filepath = 'C:\Users\jsamuels\Documents\MATLAB\Thesis\Animations\';



%figure(1)
for iFrame = 1:nFrames
    %imageFileName = [ num2str(positionsInMM(iFrame)) 'mm' 'flux.mat'];
    %filename = [sliceFilePath imageFileName];
    %load(filename);  % Frame loaded to fluxMat
    
    
    % Image Processing
    pic = ((Science_Target_Frames(:,:,iFrame))./Average_Max);
    pic(pic<=0) = eps;
    fig = figure('visible', 'off'); 
    imagesc(x_cam_phys, y_cam_phys, log10(pic), [-10 -8]);
    colorbar; axis square;
    %xlim([-OPA OPA]); ylim([-OPA OPA]);
    xlabel('Millimeters');
    title('CCD Expectation(physical units)');
    set(gca, 'Ydir', 'normal');
    savefig('CCD_Expectation');
    
    
    %Labeling
    
    
    M(iFrame) = getframe(fig);
    close all;
end

filename = [filepath videoName '.avi'];
%filename = [videoName '.avi'];
%movie2avi(M, filename,'fps', 4);
% Better video writer:
writerObj = VideoWriter(filename, 'Motion JPEG AVI');
writerObj.FrameRate = 2;
open(writerObj);
writeVideo(writerObj,M);
close(writerObj);
% % PCA_PSF_Subtraction simulates the PSF subtraction that would be used to
% remove unwanted noise from an optical simulation and ideally reveal
% science targets. 
%
% This function takes as input a set of reference and science target
% frames, and it then assembles them and implements PSF subtraction using
% Principal Component Analysis

% Author: Josh Samuels
% March 2015

function [Cleaned_Science_Target, Z, sigma] = PCA_PSF_Subtraction(configFile, Reference_Star_Frames, Science_Target_Frames, K_Components)
  

    load(configFile);

    % Code assumes that Ref and science frames will match in dimension,(no
    % checks)
    picHeight = size(Reference_Star_Frames,1);
    picWidth = size(Reference_Star_Frames,2);
    numRefs = size(Reference_Star_Frames,3);
    numSci = size(Science_Target_Frames,3);

    vectorizedPicLength = picHeight*picWidth;
    M = zeros(numRefs,vectorizedPicLength);

    % Vectorize input photographs into large Matrix for SVD

        % build m x n matrix of reference pics (rows) organized by pixel (columns)
    for im = 1:numRefs
        for r = 1:picHeight
            M(im,(((r-1)*picWidth)+1):(r*picWidth)) = Reference_Star_Frames(r,:,im);
        end
    end

    T = zeros(numSci,vectorizedPicLength);
        % build matrix of row Vectors T that represent science frames
    for im = 1:numSci
        for r = 1:picHeight
            T(im,(((r-1)*picWidth)+1):(r*picWidth)) = Science_Target_Frames(r,:,im);
        end
    end

    % Generate PCA Basis using Reference Star Frames

    % Subtract out mean of each column
    S_bar = mean(M,1);
    % Let R be the mean centered version of reference image matrix M
    R = M;
    for row = 1:numRefs
        M(row,:) =  M(row,:) - S_bar; 
    end

    % Compute Singular Value Decomposition:
    % U Sigma V^transpose = R
    [U,sigma,V] = svd(R);
    V_T = transpose(V);  % (V_T's rows are the "right-singular vectors of R)
    %K = 25;%floor(N/10);  % Truncate to K rows to eliminate overfitting
    K = K_Components; % Populate K given Function input
    Z = V_T(1:K,:);

    % Subtract out Noise Vectors from Science Target
    Cleaned_Science_Target = zeros(picHeight,picWidth,numSci);
    Z_T_Z = (transpose(Z))*Z;
    for im = 1:numSci
        t = T(im,:);
        % S_hat is an estimate of the PSF by expanding the science picture
        % on the basis generated from the reference frames
        S_hat = S_bar + (t-S_bar)*Z_T_Z;
        
        % Load PSF-subtracted image into output variable
        F = t-S_hat;
        for r = 1:picHeight
            Cleaned_Science_Target(r,:,im) = F(1,(((r-1)*picWidth)+1):(r*picWidth));
        end
    end
end


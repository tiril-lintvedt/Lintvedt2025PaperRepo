function [tuner] = tune_loading_dir(Mean, loading)
    % Tunes loading/eigenvector direction to correspond with original 
    % spectrum peaks direction. The correlation degree between loading and 
    % spectrum mean is checked, and provides a tuner (-1 or 1). After tuner 
    % is obtained, the tuner must be multiplied with the loading in order to 
    % flip spectrum in the direction corresponding to spectrum mean.
    % This is motivated by the inherently randomness in eigenvector directions  
    % in PCA and similar methods.
    
    % Input: spectrum mean (saisir structure), loading
    % Output: tuner (-1,1)
    
    rho1 = corrcoef(Mean.d',loading);

    % Check if upside down:
    isUpsideDown = rho1(1,2) < 0;  
 
    if isUpsideDown 
        tuner= 1 ;
    else 
        tuner = -1 ;
    end
    
    
end
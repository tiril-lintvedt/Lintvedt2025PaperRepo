function X_seg = segmentation_insets(X, insets)
% Segmentation method which discart image areas according to the pixel 
% percentages "insets"
%
% Input:
%           X       - Image in saisir format
%           insets  - List of image percentages to discart according to
%                     [top left bottom right], E.g. [0 20 20 10]
% 
% -------------------------------------------------------------------------

    img = X;
    
    % Vectorize/flatten the cube (one row correspond to one pixel): 
    cubeSize = size(single(img.d));
    mask = mean(single(img.d),3) ~= 0; % Mask to remove the background (assumed set to 0)
    
    xwidth = max(sum(mask,2));
    ywidth = max(sum(mask,1));
    
    % Insets in pixels
    inset_top = floor((insets(1)/100) * ywidth); 
    inset_left = floor((insets(2)/100) * xwidth);
    inset_bottom = floor((insets(3)/100) * ywidth);
    inset_right = floor((insets(4)/100) * xwidth);
            
    
    
    vectorizedCube = reshape(single(img.d), [cubeSize(1)*cubeSize(2), cubeSize(3)]);
    vectorizedMask = reshape(mask, [cubeSize(1)*cubeSize(2), 1]);
    
    relevantPixels = vectorizedCube(vectorizedMask == 1,:); % Non-background pixels
    
    % Transformation to Absorbance
    relevantPixelsAbsorbance = reflectance2Absorbance(relevantPixels,method); % Natural logarithm
    
    clear relevantPixels vectorizedCube mask
    
    % Represent the results of the pixel-wise processing in the original
    % spatial arrangement
    vectorizedResults = zeros(cubeSize(1)*cubeSize(2), cubeSize(3),'single'); 
    vectorizedResults(vectorizedMask==1,:) = relevantPixelsAbsorbance;
    
    clear vectorizedMask relevantPixelsAbsorbance
    
    cubeAbsorbance = reshape(vectorizedResults, [cubeSize(1),cubeSize(2),cubeSize(3)]);  
    img.d = cubeAbsorbance;
    
    clear cubeAbsorbance


end
function img = pixel_reflectance2abs(X, method)
% Transforms the reflectance images in the cell array X (containing one 
% structure per image), pixel by pixel, to absorbance.


    img = X;
    
    % Vectorize/flatten the cube (one row correspond to one pixel): 
    cubeSize = size(single(img.d));
    mask = mean(single(img.d),3) ~= 0; % Mask to remove the background (assumed set to 0)
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
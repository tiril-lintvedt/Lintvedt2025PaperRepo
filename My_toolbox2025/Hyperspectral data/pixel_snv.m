function   X_prep = pixel_snv(X)
% Given a cell array of structures (one per image), this function preprocesses
% each non-background pixel by SNV. Background pixels are assumed
% to be set to 0. Output is a new cell array of structures (one per image)

X_prep = {};

for i_image = 1:length(X)
    
    img = X{i_image};
    cube = img.d ;
    name = img.i;
    
    % Vectorize/flatten the cube (one row correspond to one pixel): 
    cubeSize = size(cube);
    mask = mean(cube,3) ~= 0; % Mask to remove the background (assumed set to 0)
    vectorizedCube = reshape(cube, [cubeSize(1)*cubeSize(2), cubeSize(3)]);
    vectorizedMask = reshape(mask, [cubeSize(1)*cubeSize(2), 1]);

    relevantPixels = vectorizedCube(vectorizedMask == 1,:); % Non-background pixels
    
    % SNV on non-background pixels
    relevantPixelsSnv = snv(relevantPixels);
    
    % Represent the results of the pixel-wise processing in the original
    % spatial arrangement
    vectorizedResults = zeros(cubeSize(1)*cubeSize(2), cubeSize(3)); 
    vectorizedResults(vectorizedMask==1,:) = relevantPixelsSnv;
    cubeSnv = reshape(vectorizedResults, [cubeSize(1),cubeSize(2),cubeSize(3)]);
    
    imgSnv.d = cubeSnv;
    imgSnv.v = img.v;
    imgSnv.i = name;
    
    X_prep{i_image} = imgSnv;
    
end




end
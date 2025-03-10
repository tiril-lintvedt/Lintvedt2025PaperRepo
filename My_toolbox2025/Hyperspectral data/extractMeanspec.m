function X_meanspec = extractMeanspec(X)
% Given a cell array of structures (one per image), this function calculates 
% the mean across all non-bbackground pixel. Background pixels are assumed
% to be set to 0. Output is a structure Xmeanspec with data field containing
% a 2D matrix (one mean spectrum per sample)

% Edit: mean function now ignores nan-values

X_meanspec = {};

for i_image = 1:length(X)
    
    img = X{i_image};
    cube = img.d ;
    name = img.i;
    
    % Data vectorization for pixel-wise processing
    
    % Vectorize/flatten the cube (one row correspond to one pixel): 
    cubeSize = size(cube);
    mask = mean(cube,3) ~= 0; % Mask to remove the background (assumed set to 0)
    vectorizedCube = reshape(cube, [cubeSize(1)*cubeSize(2), cubeSize(3)]);
    vectorizedMask = reshape(mask, [cubeSize(1)*cubeSize(2), 1]);

    relevantPixels = vectorizedCube(vectorizedMask == 1,:);
    
    % Use non-background pixels to calculate a mean spectrum   
    mean_spectrum = mean(relevantPixels, 'omitnan');

    X_meanspec.d(i_image,:) = mean_spectrum;
    X_meanspec.i(i_image,:) = name;
    X_meanspec.v = img.v;   
    
end



end 
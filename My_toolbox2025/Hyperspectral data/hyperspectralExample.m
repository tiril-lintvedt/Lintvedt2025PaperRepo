clear all;
% TODO: Here you can prepare some code to iterate over all the files. I
% usually use the function dir to do that. However, if you need help with
% that do not hesitate to ask:

%% 1) Load the data and the wavelengths:
fileName = 'C:\tmp\Tiril_MATLAB\Run1_Salmon-01 v2\Run1_Salmon-01\Run1_Salmon-01.mat';
load(fileName)
fileName = 'C:\tmp\Tiril_MATLAB\Run1_Salmon-01 v2\Run1_Salmon-01\Run1_Salmon-01_layers.mat';
load(fileName)
Wavelength = str2num(Wavelength);

% Since each datafile exported from Breeze has a different name, I just
% save the contents of the cube in a variable called cube to allow the code
% to be the same for differnt data files:
cube = Run1_Salmon_01;

% Simple data visualization
% Visualize one band:
imagesc(cube(:,:,50));
% Visualize the spectra of a single pixel:
plot(Wavelength, squeeze(cube(82,82,:)),'red');

% If you are interested in a smaller region of the spectra, you can use the
% following code:
minBand = 1600;
maxBand = 1800;

minBandIndex = find(Wavelength > minBand);
minBandIndex = minBandIndex(1);
maxBandIndex = find(Wavelength > maxBand);
maxBandIndex = maxBandIndex(1);

cube = cube(:,:,minBandIndex:maxBandIndex);
Wavelength = Wavelength(minBandIndex:maxBandIndex);

% Same plot but with reduced spectral range:
plot(Wavelength, squeeze(cube(82,82,:)),'red');

%% Data vectorization for pixel-wise processing
% Vectorize the cube: 
cubeSize = size(cube);
mask = mean(cube,3) ~= 0;

vectorizedCube = reshape(cube, [cubeSize(1)*cubeSize(2), cubeSize(3)]);
vectorizedMask = reshape(mask, [cubeSize(1)*cubeSize(2), 1]);

relevantPixels = vectorizedCube(vectorizedMask == 1,:);

%% Pixel-wise processing (example of k-means clustering):
numberOfClusters = 5;
[idx,C] = kmeans(relevantPixels, numberOfClusters, 'Distance', 'cosine','MaxIter', 3000, 'Replicates',5);

%% Represent the results of the pixel-wise processing in the original
%  spatial arrangement: 
vectorizedResults = zeros(cubeSize(1)*cubeSize(2),1);
vectorizedResults(vectorizedMask==1) = idx;
spatialResult = reshape(vectorizedResults, [cubeSize(1),cubeSize(2),1]);

% Example of spatial representation of the k-means results:
imagesc(spatialResult)
imagesc(spatialResult==2)

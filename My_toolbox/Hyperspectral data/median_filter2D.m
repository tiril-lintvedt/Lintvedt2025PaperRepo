function X = median_filter2D(X_2D, kernalsize)
% Convolutional median filter for 2D arrays

% Call        :  X = median_filter2D(X  ,[3,3])
% Edge method : Ignore points close to edge according to filter radius 
% (these points are set to original value)

[nx, ny] = size(X_2D);
kx = kernalsize(1);
ky = kernalsize(2);
rx = ceil(kx/2) - 1; % radius of kernel in x direction
ry = ceil(ky/2) - 1; % radius of kernel in y direction

X = X_2D; % Set output matrix equal to X_2D to preserve values on edge points
for i = (rx+1):(nx-rx)  % Loop through non edge pixels 
   for j = (ry+1):(ny-ry)        
       X(i, j) = median(X_2D((i-rx):(i+rx), (j-ry):(j+ry)), 'all');  % Median of local patch      
   end    
end
end
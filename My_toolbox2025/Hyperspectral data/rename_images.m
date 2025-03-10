function X_renamed = rename_images(X, new_names)
% Renames the images in X according to new_names (char array) assuming 
% row-to-row correspondence between the images in X and new_names list.

X_renamed = {};

if length(X) ~= size(new_names,1)
   error('Mismatch between number of images in X and rows in new_names ') 
end

for i_image = 1:length(X)
    
    img = X{i_image};    
    IMG.d = img.d;
    IMG.v = img.v;
    IMG.i = new_names(i_image,:);
    
    X_renamed{i_image} = IMG;
end
    

end
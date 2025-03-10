function Xsel = select_images_from_identifier(X, startpos, str)
% Selects images with the given string "str" in the identifier i of each image 
% in the cell array X.

   Xsel = {};
   counter = 1;
   
   for i_img = 1:length(X)
        
        image = X{i_img};
        name = image.i;      
        
        if startpos < 0
            if startsWith(name((end+startpos):end), str)            
                Xsel{counter} = image;
                counter = counter+1;
            end     
            
        else
            
            if startsWith(name(startpos:end), str)            
                Xsel{counter} = image;
                counter = counter+1;
            end     
            
        end

   end

end
function [repetative_elements] = find_equal_cell_elements(cell_array)
% This function checks if any elements of a cell array is equal/repetative,
% and outputs a cell array where the first column contains the content of
% duplicate cell elements, and the rest of the columns contain the row numbers 
% of the original cell array, in which the respective elements are found.

repetative_elements = {};
first_equal = true;

for k = 1:numel(cell_array)
    for l = 1:numel(cell_array)
        if k~=l
            if isequal(cell_array{k},cell_array{l})            
               repetative_elements{k,1} = cell_array(l); 
               repetative_elements{k,2} = k;
               
               if first_equal == true 
                   repetative_elements{k,3} = l;
                   n=3;
               else
                   n=n+1;
                   repetative_elements{k,n} = l;                 
               end
               first_equal = false;
               
            end
        end
    end
    first_equal = true;
end



if length(repetative_elements) < 1 
   disp('No cell elements are repetative.') 
end

end
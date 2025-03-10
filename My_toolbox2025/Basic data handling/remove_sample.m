function Dataset = remove_sample(Dataset, selection)
% Removes samples defined by selection from Dataset in saisir structure. The selection can be
% either a vector/double of row numbers, a boolean list of same row size as Dataset 
% or a cell containing two elements; a string identifier followed by a 
% corresponding string position in the id tag {'selection_id',pos1:pos2}.
%
%   Input: Dataset (saisir struct)
%          selection (vector of doubles or cell)
%
%   Output: Dataset (saisir struct)


if isa(selection,'double')
    Dataset.d(selection,:) = [];
    Dataset.i(selection,:) = '';
    Dataset.v = Dataset.v;
end

if isa(selection,'cell')
    id = selection{1};
    position = selection{2};
    i = cell2mat(cellfun(@(c)any(strcmp(c(position),id)),cellstr(Dataset.i), 'UniformOutput', false));
    
    Dataset.d(i,:) = [];
    Dataset.i(i,:) = '';
    Dataset.v = Dataset.v;
    Dataset.meta(i,:) = [];    
end


end
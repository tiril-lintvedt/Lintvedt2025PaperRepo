function  gr = group_from_identifier(X,pos)
% -------- Creates group assignment from identifier ---------------------------
% Creates a vector with group assignments for samples with equal identifier
% in positions pos (For e.g use for CV segment assignment)
% -----------------------------------------------------------------------------
% Inputs:
% X       - Saisir data structure (with fields i, v and d)
% pos     - vector with positions of identifier to base groups on,[pos1:pos2]
% -----------------------------------------------------------------------------
% Output:
% gr - vector with group assignments for each sample in the saisir data 
%       e.g. [1,1,1,2,2,2,....,m]
% -----------------------------------------------------------------------------


m = size(X.i,1); % number of samples/rows in X
gr = zeros(m,1); % group assignment vector
gnum = 1;   % Group number

unique_id = char(unique(cellstr(X.i(:,pos)))); % Identify all unique id tags
nid = size(unique_id,1); % Number of unique id tags

for isample= 1:nid
    id = unique_id(isample,:);
    iequal = find_from_identifier_ind(X,pos,id); 
    gr(iequal) = gnum;
    gnum = gnum +1;
end

end
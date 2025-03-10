function X_segment_avg = average_from_identifier(X, pos)

% Averages all spectra in the saisir data structure X with same identifier
% in the given id positions "pos"'
%
%   INPUT:
%               X       -   saisir data structure
%               pos     -   vector of index positions indicating the id
%                           part 
%
% -------------------------------------------------------------------------

unique_id = char(unique(cellstr(X.i(:,pos)))); % Identify all unique id tags
nid = size(unique_id,1); % Number of unique id tags

X_segment_avg.v =  X.v ;
X_segment_avg.d =  [] ;
ID = {}; 

for isample= 1:nid
    id = unique_id(isample,:); % ID tag
    iequal = find_from_identifier_ind(X,pos,id); % Indices for samples with same id tag   
    xseg_avg = mean(X.d(iequal,:));  % Average of subset with same id tag
    
    ID(end+1) =  cellstr(X.i(iequal(1),:));
    X_segment_avg.v =  X.v ;
    X_segment_avg.d(end+1,:) =  xseg_avg;
end
whole_id =char(ID');
X_segment_avg.i = whole_id(:,pos);

end
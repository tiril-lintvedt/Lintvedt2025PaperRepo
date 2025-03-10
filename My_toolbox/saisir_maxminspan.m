function X_maxminspan = saisir_maxminspan(X)

% Compute max-min value for all spectra in the saisir data structure X with 
% same identifier
% in the given id positions "pos"'
%
%   INPUT:
%               X       -   saisir data structure
%               pos     -   vector of index positions indicating the id
%                           part 
%
% -------------------------------------------------------------------------

nsamp = size(X.d,1); % Number of unique id tags

X_maxminspan.v =  X.v ;
X_maxminspan.d =  max(X.d, [], 2) - min(X.d, [], 2) ;
X_maxminspan.i = X.i;


% for isample= 1:nid
%     X_maxminspan.d(end+1,:) =  xseg_avg;
% end


end
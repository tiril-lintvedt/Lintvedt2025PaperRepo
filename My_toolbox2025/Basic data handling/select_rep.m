function [Saisir, i_selected] = select_rep(Saisir, selection)
%  ---- Copy of the Saisir data set, ignoring the given replicate selection ----------
% Selects only one replicate out of a given number of replicates to 
% choose from. The last of the replicates is chosen. 
% The replicate spectra are defined by selection. The selection can be 
% either an integer number (2, 3, 4..) defining how many replicates there 
% are of a sample, or a vector of integers indicating the position in the 
% id string defining the replicates. In the former case, the replicate 
% samples must be placed consecutively in the sample matrix.

%
%   Input: Saisir       - Saisir struct , the data set
%          selection    - vector of doubles (integer, e.g 1:2 OR cell with identifier
%                         substring which defines replicate together with substring 
%                         position E.g. {'_2', 30:31}. , replicate
%                         selection
%
%   Output: Saisir      - Saisir struct, the new dataset, ignoring reps.
%           i_selected  - vector of row indices which were extracted from 
%                         the original data set.
% -------------------------------------------------------------------------------------

if isa(selection,'double')
    % Assumes a that there are equal number of replicates for all samples.
    startrep = selection; % Start at the last of the n replicates (might as well start at 1)
    step = selection; % number of replicates
    nsamples = size(Saisir.d,1); % Total number of aquired spectra
    
    rep_id = selection{1};
    position = selection{2};
    i = cell2mat(cellfun(@(c)any(strcmp(c(position),id)),cellstr(Dataset.i), 'UniformOutput', false));
    
    i_selected = startrep:step:nsamples; 
    Saisir.d = Saisir.d(startrep:step:nsamples,:) ;
    Saisir.i = Saisir.i(startrep:step:nsamples,:);
    Saisir.v = Saisir.v;
    warning('Assumption: equal number of replicates for all samples. ALSO: Double check output accuracy.')
   
end

if isa(selection,'cell') % Or string ?
    disp('Use function "select_from_identifier" in stead.')

    
% 
%     
%     rep_id = selection{1};
%     position = selection{2};
%     i = cell2mat(cellfun(@(c)any(strcmp(c(position),id)),cellstr(Dataset.i), 'UniformOutput', false));
%     
%     Dataset.d = Dataset.d(i,:);
%     Dataset.i = Dataset.i(i,:);
%     Dataset.v = Dataset.v;
% 
       
end



end
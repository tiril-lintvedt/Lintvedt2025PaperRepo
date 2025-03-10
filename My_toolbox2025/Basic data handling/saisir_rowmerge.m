function SaisirMerged = saisir_rowmerge(varargin)
% Input:
%           varargin - Any number of Saisir structures (to be concatenated by the data rows)

DataMerged = [];
idMerged = char; 
for i = 1:nargin
    Saisir = varargin{i};
    DataMerged = cat(1,DataMerged, Saisir.d);
    idMerged = char(vertcat(cellstr(idMerged), cellstr(Saisir.i)));
    %idMerged = vertcat(idMerged, Saisir.i);
end

SaisirMerged.d = DataMerged; %cat(1,Saisir1.d, Saisir2.d);
SaisirMerged.i = idMerged(~cellfun('isempty',cellstr(idMerged)),:) ; % Remove empty cell array due to the initiation of idMerged
SaisirMerged.v = varargin{1}.v; % Assumes that all Saisir structure have same wavenumbers  
end
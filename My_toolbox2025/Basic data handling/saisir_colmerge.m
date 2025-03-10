function SaisirMerged = saisir_colmerge(varargin)
% Input:
%           varargin - Any number of Saisir structures (to be concatenated by the data cols)

DataMerged = [];
vMerged = char; 
for i = 1:nargin
    Saisir = varargin{i};
    DataMerged = cat(2,DataMerged, Saisir.d);
    vMerged = char(vertcat(vMerged, cellstr(Saisir.v)));
end

SaisirMerged.d = DataMerged; %cat(1,Saisir1.d, Saisir2.d);
SaisirMerged.i = Saisir.i ; % Remove empty cell array due to the initiation of idMerged
SaisirMerged.v = vMerged ; % Assumes that all Saisir structure have same wavenumbers  
end
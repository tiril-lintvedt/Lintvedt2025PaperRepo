function Xsaisir = dataset2saisir(Dataset)
% ----- Converts a PLS toolbox Dataset (spectra) to Saisir structure ------
%
%   INPUT:
%               Dataset     -   A PLS Toolbox Dataset 
%               
%   OUTPUT:
%               Xsaisir     -   A Saisir data structure with fields d
%                               (data/spectra),v(frequency/xaxis) and i 
%                               (sample names). 
%
%
%   Individual data can be extracted by:
%   dataArray = Xsaisir.d;
%   frequencies = Xsaisir.v % (use function str2num and num2str to convert 
%                               between versions in char/double) 
%   sample_names = Xsaisir.i
%
%
%   NB! All other information than data values (spectra), x axis values
%       (frequency) and sample names are discarded.
%
% 
% -------------------------------------------------------------------------

Dataset = struct(Dataset);

% Exstract relevant information for spectra into Saisir structure
Xsaisir.d = Dataset.data;
Xsaisir.v = num2str(Dataset.axisscale{2,1,1}');
Xsaisir.i = Dataset.label{1,1,1};

end
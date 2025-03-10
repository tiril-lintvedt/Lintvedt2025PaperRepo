function DataTable = saisir2table(Xsaisir)
% ----- Converts a PLS toolbox Dataset (spectra) to Saisir structure ------
%
%   INPUT:
%               Dataset     -   A PLS Toolbox Dataset 
%               
%   OUTPUT:
%               DataTable     -   A Matlab table with sample names and x
%                                 axis frequencies correspondng to Saisir
%                                 data structure
%
%   NB! All other information than data values (spectra), x axis values
%       (frequency) and sample names are discarded.
% -------------------------------------------------------------------------

DataTable = array2table(Xsaisir.d,'VariableNames',string(Xsaisir.v),'RowNames',string(Xsaisir.i));

end
function [Xsaisir, Xcat] = load_unscrmat_2_saisir(fullpath)
%
% Loads a .mat file (given by full path) that was stored by unscrambler and
% transforms it into saisir format. 
%
% -------------------------------------------------------------------------


% Load -mat file
load(fullpath)

% Determine what is the data matrix variable
DoubleObj = get_objects_from_ws('double', 'caller');
if length(DoubleObj) > 1
    error('Input is not of expected Unscrambler format. More than one possible data matrix.')
else
    Data = eval(DoubleObj{1});
end


% Weed out any potential category variables
numvals = cellfun(@isempty, regexp(cellstr(VarLabels),'\w*\D\>')); % This regular expression matches matches words that do not end with a numeric digit.
othervals = ~numvals;

% Put into saisir structure
Xsaisir.i = ObjLabels;
Xsaisir.v = VarLabels(numvals,:);
Xsaisir.d = Data(:, numvals);

% Matrix of other values (likely category variables) NOW: Maybe only
% made to handle numerical category variables ?
Xcat.i = ObjLabels;
Xcat.v = VarLabels(othervals,:);
Xcat.d = Data(:, othervals);


end
function XsaisirRatios = saisir_ratio(Xsaisir, WN1, WN2)
% ----------------- Intensity ratio between two columns -------------------
% Calculates the intensity ratio at WN1/WN2 (x axis values) for each row in 
% the data matrix  Xsaisir.d
% -------------------------------------------------------------------------
%
%       INPUT: 
%                   Xsaisir     -   Saisir data structure
%                   WN1         -   Frequency/wavelength (in given x axis values)   
%                   WN2         -   Frequency/wavelength (in given x axis values)  
%
%       OUTPUT: 
%                   XsaisirRatios   -   Saisir data structure containing
%                                       the given ratio for each row corre-
%                                       sponding to Xsaisir.
% -------------------------------------------------------------------------

XsaisirRatios = Xsaisir;  % Copy for corresponing sample names

[~,i1] = saisir_getxindex(Xsaisir, WN1);
[~,i2] = saisir_getxindex(Xsaisir, WN2);

XsaisirRatios.v = ['Ratio :', num2str(WN1),'/', num2str(WN2)];
XsaisirRatios.d = Xsaisir.d(:,i1)./Xsaisir.d(:,i2);


end
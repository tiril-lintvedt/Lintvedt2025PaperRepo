function [saisir]=saisir_derivative(saisir1,K,F,Dn)
% ----- N-th order derivative using the Savitzky-Golay coefficients -------
% [saisir]=saisir_derivative(saisir1,polynome_order,window_size,derivative_order)
%
%   INPUT: 
%           saisir1      -      Spectra, SAISIR format, A column vector, or  
%                               matrix with nchannels*nspectra matrix
%           K            -      Polynomial order    
%           F            -      Window size
%           Dn           -      Derivative order
%
%   OUTPUT: 
%           Y_c          -      Saisir data structure, the corrected spectra 
%                               (baseline removed / spectra smoothed)  
%           baseline     -      Matrix, baselines for each spectrum
%           wgts         -      Matrix, weights for each spectrum
%
% example: res=saisir_derivative(x,3,21,2);
% Compute the second derivative using a polynom of power 3 as model
% and a window size of 21
%
% Reference: 
% C.B.Y. Cordella, D. Bertrand, SAISIR: A new general chemometric toolbox, 
% Trends in Analytical Chemistry 54 (2014) 75–82. http://dx.doi.org/10.1016/j.trac.2013.10.009. 
%
% Savitzky, A., Golay, M.J.E., 1964. Smoothing and Differentiation of Data 
% by Simplified Least Squares Procedures. Anal. Chem. 36, 1627–1639. https://doi.org/10.1021/ac60214a047
% -------------------------------------------------------------------------


[~,g]=sgolaycoef(K,F);
[nrow,ncol]=size(saisir1.d);
M=zeros(nrow,ncol);

F1 = (F+1)/2;
F2 = -F1+1:F1-1;

for j = F1:ncol-(F-1)/2 %Calculate the n-th derivative of the i-th spectrum
    if Dn == 0
        z = saisir1.d(:,j + F2)*g(:,1);
    else
        z = saisir1.d(:,j + F2)*(Dn*g(:,Dn+1));
    end
    M(:,j) = z;
end

saisir.d=M;
saisir.i=saisir1.i;
saisir.v=saisir1.v;





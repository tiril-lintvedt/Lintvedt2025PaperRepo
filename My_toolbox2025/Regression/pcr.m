function [b0, B] = pcr(X,y, A)

% ------- implementation of the Principal Component Analysis of data-matrix X -------
% Inputs:
% X - data matrix.
% A - the maximal number of PCs to be extracted.
% y - corresponding response vector
% -------------------------------------------------------------------------
% The function returns:
% b0 - the constant terms (1×A vector) of the PCR-models for up to A components.
% B  - the PCR regression coeffs (m×A matrix) for the X-data for up to A components.
% -------------------------------------------------------------------------
% Additional explanations:
% V  - the PCA-loadings (conventionally called P in other litterature)
% U  - the normalized principal components
% sigma - singular values from the SVD
% -------------------------------------------------------------------------
% Predictions on new observation x: yhat = b0 + x*B
% -------------------------------------------------------------------------
% Reference: MATH310 course (Ulf Indahl) 
% -------------------------------------------------------------------------
[m,n] = size(X);
A  = min(A, min(n,m)-1);  % Make sure that the maximum number of components is not too large.

mx  = mean(X,1); % X column mean values
my = mean(y,1); % Response vector mean 
X = X - mx; % centered X data
y = y - my; % centered y data

[U, sigma, V] = svd(X,'econ'); % Singular Value Decomposition of Centered X (Thin version)
sigma = diag(sigma); % Get diagonal elements
U = U(:,1:A); sigma = sigma(1:A); V = V(:,1:A); % Restriction to the first mc components.
q = sigma.^(-1).*(U'*y);   % The regression coeffs for the PCA-scores.
B = cumsum(V.*q',2); % the PCR-regression coeffs for the X-data based on up to max A components.
b0 = my - mx*B;      % the correspondning constant terms for the PCR-models.


end

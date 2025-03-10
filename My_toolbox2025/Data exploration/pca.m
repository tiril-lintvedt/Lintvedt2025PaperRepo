function [T,cum_var,var,mx,V,U] = pca(X, A)

% ------- implementation of the Principal Component Analysis of data-matrix X -------
% Inputs:
% X  - data-matrix.
% A - the maximal number of PCs to be extracted.
% -------------------------------------------------------------------
% The function returns:
% T  - matrix of PCA-scores
% var - the PCA-variances (in percent of total variance)
% cum_var - The cumulative PCA variances (In percent of total variance)
% mx - row vector of the X-column mean values
% V  - the PCA-loadings (conventionally called P in other litterature)
% U  - the normalized principal components
% -------------------------------------------------------------------------
% Reference: MATH310 course (Ulf Indahl)
% -------------------------------------------------------------------------
[m,n] = size(X);
A  = min(A, min(n,m)-1);  % Make sure that the maximum number of components is not too large.
mx  = mean(X,1);
[U, sigma, V] = svd(X - mx,'econ'); % Singular Value Decomposition of Centered X (Thin version)
sigma = diag(sigma); % Get diagonal elements
var_tot = sum((sigma.^2)./m); % total variance
U = U(:,1:A); sigma = sigma(1:A); V = V(:,1:A); % Restriction to the first mc components.
T   = U.*sigma' ;     % the PCA-scores
var  = (sigma.^2)./m ;  % the PCA-variances
var = (var./ var_tot).*100; % the PCA-variances in percent of total variance
cum_var = cumsum(var); % the cumulative PCA-variances expressed in percent of total variance

end

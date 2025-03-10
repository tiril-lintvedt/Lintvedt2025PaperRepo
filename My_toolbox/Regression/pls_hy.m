function [b0,B,T,W,P,q] = plsHy(X,y,A)
% ------------- My implementation of the hybrid PLS with reorthogonalization and y-deflation -------------------
% Inputs:
% X  - data-matrix.
% y  - corresponding response vector.
% A - the maximal number of PLS components to be extracted.
% -------------------------------------------------------------------
% The function returns:
% b0 - the constant terms (1×A vector) of the PLS-models for up to mc components.
% B  - the PLS regression coeffs (m×mc matrix) for the X-data for up to mc components.
% T  - matrix of PLS-scores.
% W  - matrix of PLS-weights.
% P  - matrix of PLS-loadings.
% q  - vector of regression coeffs for the PLS scores.
% -----------------------------------------------------------------------------
% Reference: 'Fast and stable partial least squares modelling: A benchmark 
%             study with theoretical comments' 
% (https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/cem.2898)
% -----------------------------------------------------------------------------
[m,n] = size(X);
A  = min(A, min(n,m)-1); % Assure that the number of extracted components is consistent with the size of the problem.
T   = zeros(m,A); 
W   = zeros(n,A); 
P   = zeros(n,A);
q   = zeros(1,A);  % - the regression coeffs for the PLS-scores.
my = mean(y,1); % - the mean of the response values y
mx = mean(X,1); % - row vector of the X-column mean values
X = X - mx; % - centered X data
y = y - my; % - centered response vector

B = zeros(A,2); % B stored by diagonals
beta = zeros(n,A); % Regression coefficients
w0 = X'*y ; w = w0/norm(w0); 
t = X*w; rho= norm(t); t = t/rho; q(1) = y'*t;
W(:,1) = w; T(:,1) = t; B(1,1) = rho; P(:,1) = X'*t;
d = w/rho; beta(:,1) = (t'*y)*d;
y = y - t*q(1);
% ---------------- Solution of the PLS1-problem ----------------
for a= 2:A
    w1 = X'*y; w = (w0-w1)/q(a-1) - rho*w; w0 = w1;  % w =X'*t - rho*w;    
    w = w - W*(W'*w); theta = norm(w); w = w/theta;  % Reorthogonalize and normalize w   
    t = X*w; t = t - T*(T'*t);  % Reorthogonalize t  % (t = X*w - theta*t)
    rho = norm(t); t = t/rho; q(a) = y'*t;
    W(:,a) = w; T(:,a) = t;
    B(a-1,2) = theta; B(a,1) = rho;
    P(:,a) = X'*t;
% ---------------- Update regression coefficients ----------------
    d = (w - theta*d)/rho ; 
    beta(:,a) = beta(:,a-1) + q(a)*d;
    y = y -t*q(a);                  % Deflate y
    
end
b0 = my - mx*beta; % The corresponding constant terms for the PLS models
B = beta; % Rename according to usual convention for PLS algorithm
end
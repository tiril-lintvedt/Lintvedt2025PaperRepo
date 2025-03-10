function [b0,B,T,W,P,q] = bidiag2(X,y,A)
% ------------- My implementation of the PLS BIDIAG2 algorithm -------------------
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
w = X'*y ; w = w/norm(w); W(:,1) = w;
t = X*w; rho= norm(t); t = t/rho; T(:,1) = t;
q(1) = y'*t;
B(1,1) = rho;
P(:,1) = X'*t;
d = w/rho; beta(:,1) = (t'*y)*d;

% ---------------- Continue bidiagonalization ----------------
for a = 2:A
    w = X'*t - rho*w; w = w - W*(W'*w);  % Reorthogonalize w 
    theta = norm(w); w = w/theta; W(:,a) = w;
    t = X*w - theta*t; t = t - T*(T'*t);  % Reorthogonalize t
    rho = norm(t); t = t/rho; T(:,a) = t;
    yt = y'*t ; q(a) = yt(1);
    B(a-1,2) = theta; B(a,1) = rho;
    P(:,a) = X'*t;
% ---------------- Update regression coefficients ----------------
    d = (w - theta*d)/rho;
    beta(:,a) = beta(:,a-1) + (t'*y)*d;
end
b0 = my - mx*beta;  % The corresponding constant terms for the PLS models
B = beta;  % Rename according to usual convention for PLS algorithm 
end
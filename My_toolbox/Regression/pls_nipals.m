function [b0,B,T,U,W,P,q] = pls_nipals(x,y,A)
% ------------- My implementation of the PLS NIPALS algorithm -------------------
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

[m,n] = size(x);
A  = min(A, min(n,m)-1); % Assure that the number of extracted components is consistent with the size of the problem.
T   = zeros(m,A); % X scores 
U   = zeros(m,A); % Y scores, NB! For Y with only 1 col, this is simply centered y values ( with deflation for the given comp)
W   = zeros(n,A); 
P   = zeros(n,A);
q   = zeros(1,A);    % - the regression coeffs for the PLS-scores.

my = mean(y,1); % Response vector mean
mx = mean(x,1); % X column mean values
X = x - mx; % centered X data
y = y - my; % centered y data

for a = 1:A
    w = X'*y;  w = w/norm(w);  W(:,a) = w;      % X'*y ~ covariances between X cols and y. 
    t = X*w;  t = t/norm(t);  T(:,a) = t; 
    yt = y'*t ; qnorm = yt /norm(yt) ;  % double check
    u = y*qnorm; U(:,a) = u;            % double check
    
% ------------------- Deflate X and y ----------------------
    P(:,a) = X'*t;         X = X - t*P(:,a)';
    yt = y'*t; q(a)   = yt(1);    y = y - q(a).*t;
end
% ---------- Calculate regression coefficients -------------
B = cumsum(W*(pinv(P'*W)*diag(q')),2); % Regression coefficient
b0 = my - mx*B; % Constant term for the regression model
end
function [k,R]=MLR(X,Y)
% ---------------------- Multiple linear regression -----------------------
% 
%   INPUT:
%
%   OURPUT:
%                k  - Coefficients of the linear regression found by 
%                       solving the linear system y = V*k
%
%                R  - Regression residuals
% 
% -------------------------------------------------------------------------

V=[ones(size(X,1),1), X]; % Vandermonde matrix
k=(V'*V)\V'*Y;  % Coefficients k of the linear regression found by solving the linear system y = V*k

YP=X*k(2:end)+k(1);
R=Y-YP;

end
function [sr,srFcrit] = SR(X,beta)
%% Selectivity Ratio
% sr = SR(X,beta);
% Assume centered data X

Ttp = X*(beta./norm(beta));
Xtp = Ttp*((X'*Ttp)/(Ttp'*Ttp))';
Xr  = X-Xtp;
sr  = sum(Xtp.*Xtp)./sum(Xr.*Xr);
sr(isnan(sr)) = 0;

% Critical threshold value found by assessing the  SR against  F
% distribution with n-2 and n-3 degrees of freedom (suggested by Rajalahti et al, 2009)
alpha = 0.05;
n = size(X,1);
srFcrit = ftest(alpha,n-2,n-3); 

% Critical threshold: sr > srFcrit

end



function bias = cal_bias(y, ypred)
% --------------- Calculation of bias of a linear model ------------------
% Reference: Bellon-Maurel et al. 2010, https://doi.org/10.1016/j.trac.2010.05.006

bias = mean(ypred) - mean(y);

end
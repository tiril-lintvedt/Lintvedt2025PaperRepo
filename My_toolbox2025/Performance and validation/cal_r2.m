function R2 = cal_r2(y, ypred)
% ---------------- Coefficient of Determination --------------------------------
% Calculates the R squared metric for a model which yields predictions
% ypred, given the real values y. Ypred can contain multiple coumns of
% predicted values (i.e from PLSR models with many different number of 
% components included)

ymean = mean(y);
SS_tot = sum((y - ymean).^2); 
SS_res = sum((y - ypred).^2);

R2 = 1 - SS_res/SS_tot ;

end
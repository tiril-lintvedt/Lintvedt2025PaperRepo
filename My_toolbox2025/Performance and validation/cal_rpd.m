function  RPD = cal_rpd(y, ypred)
% ---------------- Ratio of Prediction to Deviation (RPD) --------------------
% Calculates the RPD performance metric of a model, given the real values 
% y and the predicted values ypred. Ypred can conatin multiple coumns of
% predicted values (i.e from PLSR models with many different number of 
% components included)

rmse = cal_rmse(y, ypred);
RPD = std(y)./rmse;

end
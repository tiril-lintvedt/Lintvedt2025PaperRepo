function rmse = cal_rmse(y,ypred)
% Calculates the Root Mean square metric of a model, given the real values 
% y and the predicted values ypred. Ypred can conatin multiple coumns of
% predicted values (i.e from PLSR models with many different number of 
% components included)

nsamples = size(y,1);
rmse = sqrt((1/nsamples)*sum((y-ypred).^2));

end
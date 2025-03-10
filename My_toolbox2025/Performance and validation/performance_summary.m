function  pmetrics = performance_summary(y, ypred, p)
% Input:
%            y          - true y values (reference)
%            ypred      - predicted y values
%            p          - number of explanatory variables in the model
% 
% Output:
%            pmetrics   - saisir struct, performance summary and y 
%                         descriptive statistics

% Performance statistics
rmse = cal_rmse(y, ypred);
r2 = cal_r2(y, ypred);
%r2adj = cal_r2adj(y, ypred, p);
r2adj = [];
rpd = cal_rpd(y, ypred);

% Decriptive statistics of y
ymin = min(y);
ymax = max(y);
ymean = mean(y);
ystd = std(y);


pmetrics.i = 'Summary of model performance metrics';
pmetrics. v = char({'RMSE';'R2';'RPD';'R2adj'; 'y min'; 'y max'; 'y mean'; 'y std'}) ;
pmetrics.d = [rmse,r2,rpd, r2adj, ymin, ymax, ymean, ystd];


end
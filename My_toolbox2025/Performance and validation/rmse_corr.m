function RMSE = rmse_corr(reference,prediction)
% -------- Slope and offset corrected Root Mean Squared Error --------

% If prediction has more than one column, separate R2 values are calculated
% for each column of predictions, using the same reference vector.
% Reference must be a one-column vector of reference values.

RMSE = []; 
N = size(reference,1); % number of samples

for i = 1:size(prediction, 2)
    [~,Ri] = MLR(prediction(:,i),reference);
    RMSE(i)= sqrt((1/N)*sum(Ri.^2));    
    
end


end
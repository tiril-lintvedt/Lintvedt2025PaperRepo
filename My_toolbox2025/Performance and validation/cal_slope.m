function ab = cal_slope(y, ypred)
% ------- Calculation of slope and intercept of linear regression --------
% y = ax + b, where a = ab(1,:) and b = ab(2,:);

ab = []; 

for i = 1:size(ypred,2)
    ab(:,i) = polyfit(y,ypred(:,i), 1)'; % Fit line (1st degree polynomial)
end

end
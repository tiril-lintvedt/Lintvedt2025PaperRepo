function residuals = residuals_corr(reference,prediction)
% ------------ Slope and offset corrected predicion residuals -------------

% If prediction has more than one column (n x m), separate residual values  
% (n x m) are calculated for each column of predictions , using the same 
% reference vector. Reference must be a one-column vector of reference values.
% 


%reference=reference(:); prediction=prediction(:);

residuals = zeros(size(prediction)); % R2
for i = 1:size(prediction, 2)
    [~,Ri] = MLR(prediction(:,i),reference);
    residuals(:,i)= Ri;    
end

end
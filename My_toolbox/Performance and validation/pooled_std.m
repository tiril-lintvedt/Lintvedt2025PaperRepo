function poolstd = pooled_std(Ypred, pos) 

%  -------------------Pooled standard deviation ---------------------------
% Average spread of all data points about their group mean (not the overall
% mean). It is a weighted average of each group's standard deviation. The 
% weighting gives larger groups a proportionally greater effect on the 
% overall estimate. Pooled standard deviations are used in 2-sample t-tests, 
% ANOVAs, control charts, and capability analysis and to calculate the 
% repeatability of predictions on replicate measurements
% -------------------------------------------------------------------------
%
%          Input:
%                   Ypred   - predicted values, saisir struct
%                   pos     - position in identification string, to define
%                             replicates from
% 
%          Output: The pooled standard deviation in the set
%
% -------------------------------------------------------------------------

samples_groups = group_from_identifier(Ypred, pos); % Sample groups
M = length(unique(samples_groups)); % Number of different samples /groups
poolstd = 0; % Initialize the pooled std

for i = 1:M
    N = sum(samples_groups == i); % Number of replicates (subsamples in group)
    group_pred = Ypred.d(samples_groups == i); % Predictions of replicates
    group_mean = mean(group_pred);
    
    for j = 1:N
        ijy = group_pred(j);
        poolstd = poolstd + ((ijy - group_mean).^2)/(M*(N-1));
    end 
end

poolstd = sqrt(poolstd);

end
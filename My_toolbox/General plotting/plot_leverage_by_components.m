function plot_leverage_by_components(X_scores, Y_scores,ncomp, names, pos)
% Plots X and Y leverages for a PLS model for all components, as aid for outlier detection
% 
% Input: 
%           Xscores     - PLS X scores (T) [nsamples x ncomp]
%           Yscores     - PLS Y scores (U) [nsamples x ncomp]
%           ncomp       - int, Upper component limit to consider
%           names       - Sample names (row-to-row correspondance w/scores)
%           pos         - int1:int2, string positions in names to be plotted

% 
% -------------------------------------------------------------------------

for icomp = 1:ncomp
    plot_leverage(X_scores, Y_scores, icomp, names, pos)
    title(['Component ', num2str(icomp)], 'FontSize', 18)
    
end
end


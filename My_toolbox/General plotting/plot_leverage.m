function plot_leverage(X_scores, Y_scores, icomp, names, pos)
% Plots X and Y leverages for a PLS model, as aid for outlier detection
% 
% Input: 
%           Xscores     - PLS X scores (T) [nsamples x ncomp]
%           Yscores     - PLS Y scores (U) [nsamples x ncomp]
%           names       - Sample names (row-to-row correspondance w/scores)
%           pos         - int1:int2, string positions in names to be plotted
%           icomp       - int, which component to plot leverages for
% 
% -------------------------------------------------------------------------

X_leverage = - mean(X_scores(:, icomp)) + X_scores(:, icomp).^2;
Y_leverage = - mean(Y_scores(:, icomp)) + Y_scores(:, icomp).^2;

scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2])
plot(X_leverage,Y_leverage, '.')
hold on
scatter_names(X_leverage, Y_leverage, names, pos)
xlim([min(X_leverage)-mean(X_leverage) max(X_leverage)+mean(X_leverage)])
ylim([min(Y_leverage)-mean(Y_leverage) max(Y_leverage)+mean(Y_leverage)])
xlabel('X leverage', 'FontSize',16)
ylabel('Y leverage','FontSize',16)
set(gcf,'Color',[1 1 1])




end
function plot_var(T,cum_var,var,mx,V,U, ncomp, wn)
% ------- Plot main results from Principal Component Analysis of data-matrix X -------
% Inputs: In accordance with function pca (from MATH310)

% T  - matrix of PCA-scores
% var - the PCA-variances (in percent of total variance)
% cum_var - The cumulative PCA variances (In percent of total variance)
% mx - row vector of the X-column mean values
% V  - the PCA-loadings (conventionally called P in other litterature)
% U  - the normalized principal components
% ncomp - Number of components to plot results for
% wn    - x axis values from saisir data structure
% -------------------------------------------------------------------------
% Reference: MATH310 course (Ulf Indahl)
% -------------------------------------------------------------------------

scrsz = get(0,'ScreenSize');
figure('Position',[880 335 scrsz(3)/2.5 scrsz(4)/2])
yyaxis left
plot(1:ncomp, cum_var(1:ncomp), '-*')
ylabel('Cumulative variance explained (%)', 'Fontsize', 14)
xlabel('Number of components', 'Fontsize', 14)
yyaxis right
plot(1:ncomp, var(1:ncomp), '-*')
ylabel('Variance explained (%)', 'Fontsize', 14)

set(gcf, 'Color', [1 1 1])



end
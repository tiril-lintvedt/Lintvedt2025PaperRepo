function plot_pca(T,cum_var,var,mx,V,U, ncomp, wn)
% ------- Plot main results from Principal Component Analysis of data-matrix X -------
% Inputs: In accordance with function pca (from MATH310)

% T  - matrix of PCA-scores
% var - the PCA-variances (in percent of total variance)
% cum_var - The cumulative PCA variances (In percent of total variance)
% mx - row vector of the X-column mean values
% V  - the PCA-loadings (conventionally called P in other litterature)
% U  - the normalized principal components
% ncomp - Number of components to plot results for
% -------------------------------------------------------------------------
% Reference: MATH310 course (Ulf Indahl)
% -------------------------------------------------------------------------

scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/2 scrsz(4)-150])

for i = 1:ncomp
    subplot(ncomp,1,i)
    plot(str2num(wn), V(:,i) )  
    xlabel ('Raman shift (cm^{-1})', 'Fontsize', 14)
    ylabel(['PC',num2str(i)], 'FontSize', 14)
    yline(0)
    l = legend(num2str(cum_var(i)), 'Location', 'northwestoutside');
    title(l,'Cum. Var explained (%)')
end
set(gcf, 'Color', [1 1 1])

% figure;
% title('Score plot')
% plot(T(:,1), T(:,2), '*')
% xlabel('Score PC1')
% ylabel('Score PC2')
% set(gcf, 'Color', [1 1 1])
% 
% figure;
% yyaxis left
% plot(1:ncomp, cum_var(1:ncomp), '-*')
% ylabel('Cumulative variance explained (%)', 'Fontsize', 14)
% xlabel('Number of components', 'Fontsize', 14)
% yyaxis right
% plot(1:ncomp, var(1:ncomp), '-*')
% ylabel('Variance explained (%)', 'Fontsize', 14)
% 
% set(gcf, 'Color', [1 1 1])



end
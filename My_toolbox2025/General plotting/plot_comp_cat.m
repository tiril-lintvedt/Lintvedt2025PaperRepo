function plot_comp_cat(T,cum_var,var,mx,V,U, ncomp, varnames)
% ------- Plot main results from Principal Component Analysis of data-matrix X -------
%
% T  - matrix of PCA-scores
% var - the PCA-variances (in percent of total variance)
% cum_var - The cumulative PCA variances (In percent of total variance)
% mx - row vector of the X-column mean values
% V  - the PCA-loadings (conventionally called P in other litterature)
% U  - the normalized principal components
% ncomp - Number of components to plot results for
% varnames  -  char array (as from saisir structure) with variable names
% -------------------------------------------------------------------------


scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/1.9 scrsz(4)-150])
varnames = categorical(cellstr(varnames)'); 

for i = 1:ncomp    
    subplot(ncomp,1,i)
    bar(varnames, V(:,i))
    ylabel(['PC',num2str(i)], 'FontSize', 14)

    l = legend(num2str(cum_var(i)), 'Location', 'northwestoutside');
    title(l,'Cum. Var explained (%)')
end
%xlabel ('Raman shift (cm^{-1})', 'Fontsize', 10)
set(gcf, 'Color', [1 1 1])

end
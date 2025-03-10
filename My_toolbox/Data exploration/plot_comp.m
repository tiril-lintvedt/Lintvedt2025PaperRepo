function plot_comp(T,cum_var,var,mx,V,U, ncomp, wn)
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
figure('Position',[50 50 scrsz(3)/1.9 scrsz(4)-150])
V_abs = abs(V);
for i = 1:ncomp
    % Find the main weighted regions
    [peaksindices,peakswn] = find_peaks_prom(V_abs(:,i)', wn,1,std(V_abs(:,i))/2,10,0);
    
    subplot(ncomp,1,i)
    plot(str2num(wn), V(:,i))
    ylabel(['PC',num2str(i)], 'FontSize', 14)
    yline(0)
    
    if ~isempty(peakswn)
        %xline(peakswn);
        text(peakswn',V(peaksindices,i),cellstr(num2str(round(peakswn',1))),'FontSize',7);
    end

    l = legend(num2str(cum_var(i)), 'Location', 'northwestoutside');
    box off
    title(l,'Cum. Var explained (%)')
end
xlabel ('Raman shift (cm^{-1})', 'Fontsize', 10)
set(gcf, 'Color', [1 1 1])

end
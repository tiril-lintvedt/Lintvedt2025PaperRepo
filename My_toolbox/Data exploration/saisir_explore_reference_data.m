function saisir_explore_reference_data(X, grouppos, idpos)

% Makes the most essential plots wrt visualization of reference data X
% including PCA, correlation plots and reference value distributions (hist.)
%
%       INPUT:
%               X           -   saisir struct, data block of refs.
%               grouppos    -   index positions in id string defining groups
%               idpos       -   index positions in id strings defining unique 
%                               sample
% -------------------------------------------------------------------------

% Add more sophisticated input parsing here with grouppos and idpos as optional 
% inputs

Kcomponents = size(X.d, 2);
if Kcomponents > 6
    Kcomponents = 7;
end

% Find overall and groupwise correlations
cors = corrcoef(X.d); % Overall correlations in full set
gr = X.i(:,grouppos);
uniquegroups = unique(cellstr(gr));
ngroups = size(uniquegroups,1);
groupcors = [];
for i = 1:ngroups
    Xgr = select_from_identifier(X,grouppos(1), uniquegroups(i,:));
    r = corrcoef(Xgr.d);
    groupcors(1:ngroups,1:ngroups,i) = r;    
end



% Reference value A vs reference value B 
[h,ax,BigAx] = gplotmatrix(X.d,[],gr,[],[],[],[],'grpbars',cellstr(X.v)');
% [MEANS,SEM,COUNTS,GNAME] = grpstats(X.d(:,1:2),gr); % In these outputs, rows = groups, cols = refernce type(e.g. fat/protein)

% Remove duplicate plots, and add correlation coefficients, overall and on
% group level
nrefs = size(X.d,2);
c = colororder;
for i = 2:nrefs
    for j = 1:(nrefs-1)
        if j >= i ;break;end
        cla(ax(i,j)) 
        text(ax(i,j),0.1,0.7,['r = ', num2str(round(cors(i,j),2))],'Units','normalized', 'FontSize',16)
        
        dy= 0;
        for k = 1:ngroups
            text(ax(i,j),0.1,0.4-dy,[uniquegroups{k},' : ', num2str(round(groupcors(i,j,k),2))],'Units','normalized', 'FontSize',10) % ,'Color', c(k,:)
            dy = dy + 0.1;
        end
    end
end


% Snippets without builtin funcs:
%set(gcf,'')
% scrsz = get(0,'ScreenSize');
% figure('Position',[50 50 scrsz(3)/1.5 scrsz(4)-150])
% hold on
% set(gcf, 'Color', [1 1 1])
% %title('Selected candidate samples', 'Fontsize', 16)
% c = 1;
% for i = 1:Kcomponents          
%     for j = 1:Kcomponents     
%         if j <= i
%             continue
%         end
%         subplot(3,ceil(Kcomponents/2),c)
%         plot(X.d(:,i),X.d(:,j), '*' )
%         text(0.05, 0.9, num2str(cormap(i,j)), 'Units', 'Normalized')
%         xlabel(X.v(i,:), 'Fontsize', 14)
%         ylabel(X.v(j,:), 'Fontsize', 14)
%         c = c+1;
%     end
% end

% % Distribution of sample compositions
% figure
% set(gcf, 'color', [1 1 1])
% hold on
% for i = 1:Kcomponents
%     subplot(2,ceil(Kcomponents/2),i)
%     histogram(X.d(:,i))
%     title(X.v(i,:))
% end

% PCA 
[T,cum_var,var,mx,V,U] = pca(X.d, Kcomponents-1);
plot_scores(T,X,grouppos,idpos)
plot_comp_cat(T,cum_var,var,mx,V,U,Kcomponents-1,X.v)

end
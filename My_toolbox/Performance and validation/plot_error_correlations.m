function plot_error_correlations(Performance,GroupID, IndicatorID, varargin)
% Function to check correlation between types of error and indicator
% variables, as followed bythe all-combinations validation.
% varargin: any number of indicator variables to check correlation with 
% (NB! Difference form in accordance with the all-combinations validation)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SETTINGS 
% -------------------------------------------------------------------------

% Main plot colors
plot_colors =[0 0.5 0.9   ; 
              0.9 0.5 0   ; 
              0.2 0.5 0.2 ; 
              0.9 0.3 0.1 ;
              0.7 0.7 0.2 ; 
              0   0   0   ;
              0.5 0.5 0.5 ;
              0.2 0.7 0.7 ;
              0.1 0.3 0.9 ;
              0.3 0.9 0.3 ;
              0.9 0.9 0.9 ;
              0.1 0.1 0.5 ;];

plot_markers  = {'^','o','square','diamond','*','x','hexagram','.',...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.'};

set(0,'defaultfigurecolor',[1 1 1])

% -------------------------------------------------------------------------
% MAIN 
% -------------------------------------------------------------------------

nIndicators = nargin-3;
figure('Position',[50 50 600 700]);
t = tiledlayout(nIndicators,3);%,'TileSpacing','compact');


for i = 1:nIndicators
    XInd = varargin{i};

    absdiffvals = 0;
    [DiffXInd_vect, DiffRMSEcorr_vect,DiffID] = check_error_correlations(Performance, XInd, GroupID, 'rmsepcorr',absdiffvals,'showPlot', 'no');
    [DiffXInd_vect, DiffBias_vect,DiffID] = check_error_correlations(Performance, XInd, GroupID, 'bias',absdiffvals,'showPlot', 'no'); %xline(0);
    [DiffXInd_vect, DiffSlope_vect,DiffID] = check_error_correlations(Performance, XInd, GroupID, 'slope',absdiffvals,'showPlot', 'no');
    
    nexttile%(i)
    hold on
    for j = 1:length(DiffRMSEcorr_vect)
        plot(DiffXInd_vect(j),DiffRMSEcorr_vect(j),plot_markers{j},'color', plot_colors(j,:)); 
        ylabel(['Diff ',IndicatorID{i}])
        xlabel('Diff RMSE_{corr}')
        
    end
    r = corrcoef(DiffXInd_vect,DiffRMSEcorr_vect);
    text(0.1, 0.1, ['r = ',num2str(round(r(1,2),2))],'units','normalized','color', 'r')
    axis(gca, 'padded')
    xlim(gca, xlim(gca) + [-1,1]*range(xlim(gca)).* 0.30)
    ylim(gca, ylim(gca) + [-1,1]*range(ylim(gca)).* 0.30)

    nexttile%(i)
    hold on
    for j = 1:length(DiffBias_vect)
        plot(DiffXInd_vect(j),DiffBias_vect(j),plot_markers{j},'color', plot_colors(j,:)); 
        ylabel(['Diff ',IndicatorID{i}])
        xlabel('Diff Bias')
        %xline(0)
    end
    r = corrcoef(DiffXInd_vect,DiffBias_vect);
    text(0.1, 0.1, ['r = ',num2str(round(r(1,2),2))],'units','normalized','color', 'r')
    axis(gca, 'padded')
    xlim(gca, xlim(gca) + [-1,1]*range(xlim(gca)).* 0.30)
    ylim(gca, ylim(gca) + [-1,1]*range(ylim(gca)).* 0.30)

    nexttile%(i)
    hold on
    for j = 1:length(DiffSlope_vect)
        plot(DiffXInd_vect(j),DiffSlope_vect(j),plot_markers{j},'color', plot_colors(j,:)); 
        ylabel(['Diff ',IndicatorID{i}])
        xlabel('Diff Slope')
        %xline(0)
    end
    r = corrcoef(DiffXInd_vect,DiffSlope_vect);
    text(0.1, 0.1, ['r = ',num2str(round(r(1,2),2))],'units','normalized','color', 'r')
    axis(gca, 'padded')
    xlim(gca, xlim(gca) + [-1,1]*range(xlim(gca)).* 0.30)
    ylim(gca, ylim(gca) + [-1,1]*range(ylim(gca)).* 0.30)

end

legend(GroupID,'Location','southoutside')

end
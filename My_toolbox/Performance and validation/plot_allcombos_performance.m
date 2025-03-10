function [FIG, axeshandles] = plot_allcombos_performance(Performance, varargin)
% Plots main performance metrics from the Performance struct obtained by
% the all combinations validations

% -------------------------------------------------------------------------
% SETTINGS 
% -------------------------------------------------------------------------

plot_colors =[0.1 0.1 0.5 ; 
              0.5 0.5 0.5 ;
              0.9 0.5 0   ; 
              0.2 0.7 0.7 ;
              0.6 0.2 0.9 ;
              0.3 0.8 0.3 ;
              0.9 0.3 0.6 ; 
              0   0   0   ;
              0 0.5 0.9   ; 
              0.2 0.5 0.2 ; 
              0.7 0.7 0.2 ; 
              0.9 0.3 0.1 ;
              0.1 0.1 0.5 ; 
              0.5 0.5 0.5 ;
              0.9 0.5 0   ; 
              0.2 0.7 0.7 ;
              0.6 0.2 0.9 ;
              0.3 0.8 0.3 ;
              0.9 0.3 0.6 ; 
              0   0   0   ;
              0 0.5 0.9   ; 
              0.2 0.5 0.2 ; 
              0.7 0.7 0.2 ; 
              0.9 0.3 0.1 ;
];

plot_markers  = {'^','o','square','diamond','*','x','hexagram','.',...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.',...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.'};

set(0,'defaultfigurecolor',[1 1 1])

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultMetric = {'rmsepcorr','bias','slope'}; 
defaultRmsepLim = [];
defaultBiasLim = [];
defaultSlopeLim = [];

p = inputParser;
   expectedMetrics = {'rmsep','rmsepcorr','bias','slope'};
   
   addRequired(p,'Performance'); % positional arg
   addParameter(p,'Metric',defaultMetric); % Name Value pair
   addParameter(p,'RmsepcorrLim',defaultRmsepLim); % Name Value pair
   addParameter(p,'BiasLim',defaultBiasLim); % Name Value pair
   addParameter(p,'SlopeLim',defaultSlopeLim); % Name Value pair

   parse(p,Performance, varargin{:});
   
Metric = p.Results.Metric;
% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

scrsz = get(0,'ScreenSize');
FIG = figure('Position',[50 50 scrsz(3)/2.5 scrsz(4)/1.9]);
nmetrics = length(Metric);
axeshandles = {};

tiledlayout(nmetrics,1,'TileSpacing','Compact');

for m = 1:length(Metric)
    metric = Metric{m};
    nexttile
    grid on
    hold on 
    ci = 1; 
    for i = 1:size(Performance.(metric),1) 
        for j = 1:(size(Performance.(metric),2))
            if i == j; continue;end
            %descr = [Performance.i(i,1:end-15) Performance.v(j,1:end-12)];
            descr = [char(strrep(cellstr(Performance.i(i,:)),'->','')) Performance.v(j,:)];
            b = bar(categorical(cellstr(descr)), Performance.(metric)(i,j),'facecolor', plot_colors(ci,:),'EdgeColor', plot_colors(ci,:));            
            xtips2 = b(1).XEndPoints;
            labels2 = string(round(b(1).YData,2));
            text(xtips2,Performance.(metric)(i,j)/1.9,labels2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.95 0.95 0.95],'FontSize',9)
            if strcmp(metric,'rmsepcorr')   
                text(xtips2,Performance.(metric)(i,j)+0.3,['LV=',num2str(Performance.lv(i,j))],'HorizontalAlignment','center','VerticalAlignment','top','Color',[0 0 0],'FontSize',8)
                plot(xtips2,Performance.(metric)(i,j)+0.4, 'color',plot_colors(ci,:),'Marker',plot_markers{ci})
            end

            if strcmp(metric,'slope')
                yline(1,'--')
            end
            ci = ci +1;

        end        
    end
    ylabel(metric)
     
    if m ~= nmetrics
        set(gca, 'XTickLabel',{''})
    end
    
    axeshandles{end+1} = gca;

end





end
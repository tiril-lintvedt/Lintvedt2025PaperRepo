function [DiffScattind_vect, Diffperformance_vect, description] = check_error_correlations(Performance, X_sc_avg, names, metric, absolute_difference, varargin)
% example input
% names = ['fin0', 'grov', 'hom'];

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
% INPUT PARSING
% -------------------------------------------------------------------------

defaultShowPlot = 'yes';

p = inputParser;
   expectedShowPlot = {'yes', 'no'};
  
   addRequired(p,'Performance'); % positional arg
   addRequired(p,'X_sc_avg'); % positional arg
   addRequired(p,'names'); % positional arg
   addRequired(p,'metric'); % positional arg
   addRequired(p,'absolute_difference'); % positional arg
   addParameter(p,'showPlot',defaultShowPlot,  @(x) any(validatestring(x,expectedShowPlot))); % Name Value pair
   
   parse(p,Performance, X_sc_avg, names, metric, absolute_difference, varargin{:});
   show_plot = p.Results.showPlot;


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

Metric = ['Performance.',metric];
DiffScattind_vect  = [];
Diffperformance_vect = [];
description = {};

if strcmp(show_plot, 'yes')
    figure; hold on
end

c = linspecer(size(Performance.rmsep,1)*size(Performance.rmsep,2));
phandles = {};
ci = 1;
for i =1:3
    for j = 1:3
        if i == j; continue;end;
        DiffScattind(i,j) = X_sc_avg.d(i) - X_sc_avg.d(j) ;
        
        if absolute_difference
            DiffScattind_vect(end+1) = abs(X_sc_avg.d(i) - X_sc_avg.d(j));
            Diffperformance_vect(end+1) = abs(eval([Metric,'(i,j)']));
            if strcmp(show_plot, 'yes')
                phandles{ci} = plot(abs(eval([Metric,'(i,j)'])),abs(DiffScattind(i,j)),plot_markers{ci},'color', plot_colors(ci,:));
            end
        else
            DiffScattind_vect(end+1) = X_sc_avg.d(i) - X_sc_avg.d(j);
            Diffperformance_vect(end+1) = eval([Metric,'(i,j)']);
            if strcmp(show_plot, 'yes')
                phandles{ci} = plot(eval([Metric,'(i,j)']),DiffScattind(i,j),plot_markers{ci},'color', plot_colors(ci,:));
            end
        end

        
        %text(eval([Metric,'(i,j)']),DiffScattind(i,j),[Performance.i(i,:),Performance.v(j,:)])
        description{ci} = [Performance.i(i,:),Performance.v(j,:)];
        
        if strcmp(show_plot, 'yes')
            xlabel(metric,'FontSize', 18)
            yline(0)
        end
        ci = ci+1;
    end
end

% Check correlation
r = corrcoef(DiffScattind_vect, Diffperformance_vect);

if strcmp(show_plot, 'yes')
    ylabel('Difference in indicator variable','FontSize',18) 
    text(0.1, 0.1, ['r = ', num2str(round(r(1,2),2))],'Unit', 'normalized', 'Fontsize', 14,'color', 'r')
    legend([phandles{:}],description{:})
end


end
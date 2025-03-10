function plot_scores(T,Xsaisir, colour_pos, name_pos, varargin)

% --------- Score plot for PCA --------------------------------------------
% Inputs:
%           T            - Matrix of PLS-scores.
%           colour_pos   - Position of id tag to determine group colouring
%           Xsaisir      - The original X data block (saisir data struct)
%           nametags     - Position in name tag (Whether to show nametags
%                           in plot)
%                           
%
% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultPlotNames = 'no'; 
defaultcolorBy = nan;
defautNamePos = nan;

p = inputParser;
   validLen = @(x)  (length(x)==size(T,1));
   expectedPlotNames = {'yes','no'};
   
   addRequired(p,'T'); % positional arg
   addRequired(p,'Xsaisir'); % positional arg
   addRequired(p,'colour_pos'); % positional arg
   addOptional(p,'name_pos',defautNamePos);
   addParameter(p,'PlotNames',defaultPlotNames,@(x) any(validatestring(x,expectedPlotNames)));
   addParameter(p,'colorBy',defaultcolorBy,validLen);
     
   parse(p, T, Xsaisir, colour_pos, name_pos, varargin{:});
   
   plot_names = p.Results.PlotNames; 
   colorByVec = p.Results.colorBy; 

% -------------------------------------------------------------------------


unique_id = char(unique(cellstr(Xsaisir.i(:,colour_pos)))); % Identify all unique id tags in position colour_pos
colour_groups = group_from_identifier(Xsaisir, colour_pos);
unique_groups = unique(colour_groups);
colour = linspecer(length(unique_groups), 'qualitative');
ncomps = size(T,2);

% Plot technicals
nsubp = nchoosek(ncomps,2); % Total number of subplots
nrows = ceil(nsubp/3); % Number of columns for subplots
ncols = ceil(nsubp/nrows); % Number of rows for subplots
phandles = []; 
iplot = 1;

scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)-150 scrsz(4)-150])
title('Score plot')

for s1 = 1:ncomps
    for s2 = 1:ncomps
        
        if s2 <= s1   % avoid redundant plots
           continue 
        end
        
        subplot(nrows,ncols,iplot); iplot = iplot +1;
        hold on

        if ~isnan(colorByVec)
            for i_s = 1:size(T,1)
                plot(T(i_s,s1), T(i_s,s2), '*');
            end
            color_by(colorByVec)

        else
            for i = 1:length(unique_groups)
                gr_id = find(colour_groups == unique_groups(i));
                p_i = plot(T(gr_id,s1), T(gr_id,s2), '*', 'Color',colour(i,:));
                phandles(i,1) = p_i(1);
                
                if strcmp(plot_names,'yes')
                   scatter_names(T(gr_id,s1), T(gr_id,s2), Xsaisir.i(gr_id,:), name_pos) 
                end
                
            end
        end
        
        xline(0)
        yline(0)
        xlabel(['Score PC',num2str(s1)], 'Fontsize', 14)
        ylabel(['Score PC',num2str(s2)], 'Fontsize', 14)

        % if length(unique_groups) ~= size(Xsaisir.d,1)
        %     legend(phandles,cellstr(num2str(unique_id)))
        % end
    end
    
end


if isnan(colorByVec)
    leg = legend(phandles,cellstr(num2str(unique_id)));
    title(leg, 'Group', 'Fontsize', 14)
end
set(gcf, 'Color', [1 1 1])



end


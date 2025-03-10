function phandles = plot_from_identifier2(Xsaisir,id_pos,varargin)
% ----- Plots X (row-wise) in different colors depending on identifier -----
%       (data rows with same identifier get same color)
%
%   INPUT:
%              Xsaisir   -   Saisir data structure with data point
%               id_pos   -   Positions (cols) in identifier strings to use 
%             plottype   -   The usual Matlab color/line style quick-codes
%
%
% -------------------------------------------------------------------------
% SETUP
% -------------------------------------------------------------------------

bar_colors = ['r'; 'b'; 'g'; 'k'; 'w' ;'m'; 'y'; 'c';];

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultPlotType = 'spectrum';
defaultLineSpec = '-';

p = inputParser;
   expectedPlotTypes = {'spectrum', 'idbar'};
   
   addRequired(p,'Xsaisir'); % positional arg
   addRequired(p,'id_pos'); % positional arg
   addParameter(p,'plot_type',defaultPlotType,@(x) any(validatestring(x,expectedPlotTypes)));
   addParameter(p,'line_spec',defaultLineSpec);
   parse(p,Xsaisir, id_pos, varargin{:});

   plot_type = p.Results.plot_type;
   line_spec =  p.Results.line_spec;
   
   if strcmp(plot_type,'idbar')
    col = input('Which col index num (integer) do you want to plot?');
   end

% -------------------------------------------------------------------------


% Identify the different groups based on identifier
unique_id = char(unique(cellstr(Xsaisir.i(:,id_pos)))); % Identify all unique id tags in position colour_pos
colour_groups = group_from_identifier(Xsaisir, id_pos);
unique_groups = unique(colour_groups);
ngroups = length(unique_groups);
colour = linspecer(ngroups, 'qualitative');

scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
hold on 

for i = 1:ngroups
    gr_id = find(colour_groups == unique_groups(i));

    if strcmp(plot_type, 'spectrum')
        p_i = plot(str2num(Xsaisir.v),Xsaisir.d(gr_id,:),line_spec, 'Color',colour(i,:));
    elseif strcmp(plot_type, 'idbar')
        p_i = bar(strrep(cellstr(Xsaisir.i(gr_id,:)),'_',' '), Xsaisir.d(gr_id,col), bar_colors(i,:));
        set(gcf, 'color', [1 1 1])
    end

    phandles(i,1) = p_i(1);
end

legend(phandles,unique_id, 'FontSize', 14)
set(gcf,'Color',[1 1 1])

end
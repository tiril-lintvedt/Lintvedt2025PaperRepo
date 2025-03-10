function color_by_identifier(X, idpos)
% Sets color order according to the id tag, with qualitative coloring (unique
%  id tags have same colors)
%
%           X   -   saisir data structure (with a char id field "i" )
% --------------------------------------------------------------------------

unique_id = char(unique(cellstr(X.i(:,idpos)))); % Unique ID tags
ncolours = length(unique_id);
colors = linspecer(ncolours, 'qualitative'); % Define distinguishable colours
sample_groups = group_from_identifier(X, idpos);
color_map = colors(sample_groups,:) ;
colororder(color_map)
end

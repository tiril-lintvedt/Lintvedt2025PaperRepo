function scatter_names(x,y,names, pos)
% Make a scatterplot in positions [x,y] with names.
%
%   Input:
%               x       - x positions
%               y       - y positions
%               names   - char array with names in rows
%               pos     - string position in names to plot
% -------------------------------------------------------------------------
text(x,y, strrep(cellstr(names(:,pos)),'_', ' '))

end
function cbh = color_by(x, mapName, color_range)
% Sets color order according from a value of vector representing the color
% range
% --------------------------------------------------------------------------
if ~exist('mapName','var') || isempty(mapName)
    mapName = 'jet';
end

if ~exist('crange','var') || isempty(color_range)
    color_range = [min(x) max(x)];
end

% Create color map from vector x
colour_map = vals2colormap(x, mapName); % , crange
colormap(mapName)
colororder(colour_map)

% Add color bar in plot
cbh = colorbar();
clim(color_range)


end

function phandles = plot_from_identifier(Xsaisir,id_pos,plot_type)
% ----- Plots X (row-wise) in different colors depending on identifier -----
%       (data rows with same identifier get same color)
%
%   INPUT:
%              Xsaisir   -   Saisir data structure with data point
%               id_pos   -   Positions (cols) in identifier strings to use 
%             plottype   -   The usual Matlab color/line style quick-codes
%
%
% ------------------------------------------------------------------------

if ~exist("plot_type")
    plot_type = '-';
end

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
    p_i = plot(str2num(Xsaisir.v),Xsaisir.d(gr_id,:),plot_type, 'Color',colour(i,:));
    phandles(i,1) = p_i(1);
end

legend(phandles,unique_id, 'FontSize', 14,'Location','Southeast')
set(gcf,'Color',[1 1 1])

end
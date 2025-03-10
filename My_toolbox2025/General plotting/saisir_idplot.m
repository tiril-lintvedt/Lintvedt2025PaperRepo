function hl = saisir_idplot(X,Y,I) 
% Made for spectral data, when input is only one saisir struct
% Made for scatterplot when x,y and id tag is given separately.

% X     -   Saisir data structure (or x vecor)


% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/15 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2.5])
% set(gcf,'Color',[1 1 1])
% hold on
%     
if nargin == 3 % Make scatterplot
    
    y = Y;
    x = X;
    id = strrep(cellstr(I),'_', ' ');    

    for i = 1:size(I,1)
        hl = plot(x(i),y(i),'*'); hold on
        sadf = dataTipTextRow('ID ',repelem(id(i),size(y,2)));
        hl.DataTipTemplate.DataTipRows(end+1:end+2) = sadf;
        hl.DataTipTemplate.DataTipRows(3) = hl.DataTipTemplate.DataTipRows(1);
        hl.DataTipTemplate.DataTipRows(1) = hl.DataTipTemplate.DataTipRows(4);
        hl.DataTipTemplate.DataTipRows = hl.DataTipTemplate.DataTipRows(1:3);
    end    
    xlabel('x', 'FontSize', 18)
    ylabel('y', 'FontSize', 18)
else 
    spec = X.d;
    v = str2num(X.v);
    id = strrep(cellstr(X.i),'_', ' ');

    for i = 1:size(X.i,1)
        hl = plot(v,spec(i,:),'-'); hold on
        sadf = dataTipTextRow('ID ',repelem(id(i),size(X.d,2)));
        hl.DataTipTemplate.DataTipRows(end+1:end+2) = sadf;
        hl.DataTipTemplate.DataTipRows(3) = hl.DataTipTemplate.DataTipRows(1);
        hl.DataTipTemplate.DataTipRows(1) = hl.DataTipTemplate.DataTipRows(4);
        hl.DataTipTemplate.DataTipRows = hl.DataTipTemplate.DataTipRows(1:3);
    end
    
    
end


end
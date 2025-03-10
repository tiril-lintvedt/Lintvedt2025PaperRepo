function plot_repeatability(pooled_std)

% Plots and compares the pooled standard deviation (repeatibility) across
% datasets/data blocks. 


    names = categorical(strrep(cellstr(pooled_std.i),'_', ' '));
    scrsz = get(0,'ScreenSize');
    figure('Position',[50 50 scrsz(3)/1.5 scrsz(4)/3])
    hold on
    bar(names, pooled_std.d)
    ylabel('STD_{pooled}', 'FontSize', 14)   
    set(gcf,'color', [1 1 1])


end
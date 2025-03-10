function ph = plot_Raman(wn,X)
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
    ph = plot(wn,X);
    ylabel('Intensity','FontSize',16)
    xlabel('Raman Shift (cm^{-1})', 'FontSize',16)
    set(gcf,'Color',[1 1 1])
    box off

end
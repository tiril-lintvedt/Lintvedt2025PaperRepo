function ph = plot_NIR(wl,X)
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
    ph = plot(wl,X);
    ylabel('Absorbance','FontSize',16)
    xlabel('Wavelength (nm)', 'FontSize',16)
    set(gcf,'Color',[1 1 1])
    box off

end
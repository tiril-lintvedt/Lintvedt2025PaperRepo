function plot_pixel(Zsaisir,coord, abs_trans)
% Plots spectrum at given pixel position with or without log transform

% Input:
%           Zsaisir     - Saisir image structure (with fields i, d and v)
%           coord       - Pixel coordinates, [x,y] 
%           abs_trans   - Wether or not to transform to absorbance
%                         {'none' or 'ln'}
%        

    if nargin == 2
       abs_trans = 'ln';
    end

    x = coord(1);
    y = coord(2);
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/3 50 scrsz(3)/3 scrsz(4)-150])
    if strcmp(abs_trans, 'ln')
        subplot(2,1,1)
        plot(Zsaisir.v, -log(squeeze(Zsaisir.d(y,x,:))))
        ylabel('Absorbance (ln)', 'FontSize', 14)
        xlabel('Wavelength (nm)', 'FontSize', 14)
        
        
    else
        
        subplot(2,1,1)
        plot(Zsaisir.v, squeeze(Zsaisir.d(y,x,:)))
        ylabel('Absorbance (ln)', 'FontSize', 14)
        xlabel('Wavelength (nm)', 'FontSize', 14)
        
    end
    
    subplot(2,1,2)
    imagesc(Zsaisir.d(:,:, 50))
    xlabel('x')
    ylabel('y')
    xline(x)
    yline(y)
    set(gcf, 'Color', [1 1 1])
end
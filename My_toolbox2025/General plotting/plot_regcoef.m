function FIG = plot_regcoef(wn, R,  th_scale, varargin)
%  Plotting regression vectors and detecting regions of importance
%
%       INPUT:
%                   wn          -       str, x axis values
%                   R           -       Regression vector, one or more with
%                                       one reg vector per row, with corre-
%                                       sponding x axis values in wn.
%                   th_scale    -       threshold for detection of peaks, 
%                                       given as a scaling number of the       
%                                       standard deviation of the reg.coef
%                                       i.e. detectionlimit = th_scale*STD
%                   varargin    -       Other input name/value pairs as defined
%                                       in parsing
% -------------------------------------------------------------------------
% SETTINGS
% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultThresholdScale = 0.5;
defaultScaledPLot = 'no';
defaultFindPeaks = 'yes';
defaultPlotColors = linspecer(size(R,1));

p = inputParser;
   expectedScaledPlot = {'yes','no'}; 
   expectedFindPeaks = {'yes','no'}; 

   addRequired(p,'wn'); % positional arg
   addRequired(p,'R'); % positional arg
   addOptional(p,'th_scale',defaultThresholdScale);
   addParameter(p,'scaledPlot',defaultScaledPLot,  @(x) any(validatestring(x,expectedScaledPlot)));
   addParameter(p,'findPeaks',defaultFindPeaks,  @(x) any(validatestring(x,expectedFindPeaks)));
   addParameter(p,'plotColors',defaultPlotColors, @(x) (size(x,1) >= size(R,1))); 

  
   parse(p,wn, R, th_scale, varargin{:});
   
   th_scale = p.Results.th_scale; 
   scale_plot = p.Results.scaledPlot; 
   find_peaks = p.Results.findPeaks;
   plot_colors = p.Results.plotColors;

% -------------------------------------------------------------------------
% MAIN CODE
% -------------------------------------------------------------------------

scrsz = get(0,'ScreenSize');
FIG = figure('Position',[50 50 scrsz(3)/1.5 scrsz(4)/4]);
%figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
hold on 

for i = 1:size(R,1)
    R_abs = abs(R(i,:));

    % Find the main weighted regions
    if strcmp(find_peaks,'yes')
        [peaksindices, peakswn] = find_peaks_prom(R_abs, wn,1,th_scale*std(R_abs),5,0);
    end

    if strcmp(scale_plot,'yes')
        R(i,:) = snv(R(i,:));
    end

    %subplot(ncomp,1,i)
    plot(str2num(wn), R(i,:),'Color',plot_colors(i,:))
    
    if strcmp(find_peaks,'yes') && ~isempty(peakswn)
        %xline(peakswn);
        text(peakswn',R(i,peaksindices),cellstr(num2str(round(peakswn',1))),'FontSize',7);
    end

end
yline(0) 
xlabel ('X axis values', 'Fontsize', 18)
set(gcf, 'Color', [1 1 1])
grid on
box off


if strcmp(scale_plot,'yes')
    ylabel('Reg. coefficients (SNV)', 'FontSize', 18)  
else
    ylabel('Reg. coefficients', 'FontSize', 18)   
end


end



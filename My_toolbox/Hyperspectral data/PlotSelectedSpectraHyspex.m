function raw_spectra = PlotSelectedSpectraHyspex(wl,img, varargin)
%
% function [raw_spectra, abs_spectra] = PlotSelectedSpectra(wl, img, opts)
%
% The optional structure opts can contain from none to all of the following
% fields:
%
% - range1      : integer array (1x2) 
% - range2      : integer array (1x2) 
% - autorange   : string
% - extraWidth  : integer array (1x2)
% - fig         : integer (1x1)
%
% The default values of these parameters (if not specified) are:
% range1     = [ 1 size(img,1) ]
% range2     = [ 1 size(img,2) ]
% autorange  = 'meat'
% autorange_type  = 'NIR'
% extraWidth = [ 0 0 ]
% fig        = 10
%
% If present, range1 and range2 can be used to select a range to be plotted
% along dimension 1 (vertical) and 2 (horisontal) respectively.
% 
% If the string autorange is set, it is used together with the function
% extractFood and foodFilter to autoselect a portion of the image. The
% string must be one of the known parameters of foodFilter. If options
% range1 and/or range2 is set, these override autorange.
%
% Parameter extraWidth can be used to select spectra from a wider area of
% pixels around each selected pixel. It is a vector with two elements, one
% for each dimension. The number of pixels specified is added on each side
% of the selected pixel to form a rectangle.
% 
% The parameter fig determines the figure number to use for plotting.
%
% Version Changes: 
% * 30.09.21 Tiril Lintvedt: replaced getpts function (from matlab 
%   image processing toolbox which is not available) with manual input
%   function.
%

% Handle required and optional arguments
p = inputParser;
p.addRequired('img', @isnumeric);
p.addOptional('range1', [1 size(img,1)]);
p.addOptional('range2', [1 size(img,2)], @isnumeric);
p.addOptional('autorange', 'meat', @ischar);
p.addOptional('autorange_type', 'NIR', @ischar);
p.addOptional('extrawidth', [0 0], @isnumeric);
p.addOptional('fig', 10, @isnumeric);
p.KeepUnmatched = true;
p.parse(img, varargin{:});
opts = p.Results;
explodeStruct(p.Results);

% Use autorange if neither range1 nor range2 is specified
% if numel(intersect(p.UsingDefaults, { 'range2', 'range1' })) == 2
%     bitmap  = extractFood(img, 'food', autorange, 'spectrum_type', autorange_type, opts);
%     [I,J] = find(bitmap);
%     
%     range1(1) = 1;
%     range1(2) = size(img,1);
%     
%     range1(1) = max(min(I) - 2,1);
%     range1(2) = min(max(I) + 2,size(img,1));
%     range2(1) = max(min(J) - 5,1);
%     range2(2) = min(max(J) + 5,size(img,2));
% end


% Show image
figure(fig), clf
subplot(121);
imagesc(img(range1(1):range1(2),range2(1):range2(2),20)); % 50
%imagesc(mean(img(range1(1):range1(2),range2(1):range2(2),:),3));
% imagesc(img(range1(1):range1(2),range2(1):range2(2),13)./( img(range1(1):range1(2),range2(1):range2(2),4) + 0.01 ), [0 1.5]);
title('Image','FontSize', 14)
setAxisLabels(fig,[],range1,range2);
set(fig, 'ResizeFcn', {@setAxisLabels, range1, range2});
set(fig,'Position',[0 200 1000 600])
set(gcf,'Color', [1 1 1])
grid minor

% Make the user select pixels
% Can perhaps use this to select and average spectra inlarger area...
% waitforbuttonpress
% pos = rbbox; 
% x1(1) = pos(1);x1(2)= pos(1)+ pos(3);
% x2(1) = pos(2);x2(2)= pos(2)+ pos(4);
% 
% x1 = round(x1(1):x1(end));
% x2 = round(x2(1):x2(end));
% Does not give one x and y for each point
%annotation('rectangle',pos,'Color','r') 


[x1,x2] = ginput; % Remember to push enter after selecting
x1 = round(x1);
x2 = round(x2);
labels = num2str([1:length(x1)]');

% Extract selected data from image
for i=1:length(x1)
    raw_spectra(i,:) = ExtractSpectrum(img(range1(1):range1(2),range2(1):range2(2),:),x1(i),x2(i),extrawidth);
end
text(x1,x2,labels,'FontWeight','bold','Color','white');



% Plot spectra of selected points
figure(fig), subplot(122)
plot(wl,raw_spectra, 'LineWidth', 2)
title('Spectra (not flipped)', 'FontSize',14)
legend(labels);
xlabel('Wavelength (nm)', 'FontSize',14)
ylabel('Raw signal', 'FontSize',14)
text(ones(size(x1)),raw_spectra(:,1),labels,'FontWeight','bold');



% % Plot raw spectra of selected points
% figure(fig), subplot(122)
% plot(wl,log10(1./raw_spectra)', 'LineWidth', 2)
% title('Raw spectra (not flipped)')
% legend(labels);
% text(ones(size(x1)),raw_spectra(:,1),labels,'FontWeight','bold');


% % Plot absorbance spectra of selected points
% figure(fig), subplot(223)
% abs_spectra = reflectance2Absorbance(fliplr(raw_spectra),'log');
% plot(abs_spectra', 'LineWidth', 2)
% title('Absorbance spectra')
% % legend(labels);
% text(ones(size(x1)),abs_spectra(:,1),labels,'FontWeight','bold');
% 
% % Plot SNV-normalised spectra of selected points
% figure(fig), subplot(224)
% snv_spectra = normaliseSpectra(abs_spectra, 'SNV');
% plot(snv_spectra', 'LineWidth', 2)
% title('SNV spectra')
% % legend(labels);
%text(ones(size(x1)),snv_spectra(:,1),labels,'FontWeight','bold');
%all=[raw_spectra abs_spectra snv_spectra];


function setAxisLabels(hObject, eventdata, range1, range2)
xtick = get(gca,'XTick')';
ytick = get(gca,'YTick')';
set(gca,'XTickLabel', num2str(xtick + range2(1) - 1 ));
set(gca,'YTickLabel', num2str(ytick + range1(1) - 1 ));



function spectrum = ExtractSpectrum(img,x1,x2,extrawidth)
% Utility-function for extracting a mean spectra in a given coordinate
extra1 = extrawidth(1); % Numbers of extra pixels in dimension 1
extra2 = extrawidth(2); % Numbers of extra pixels in dimension 2
start1 = max(x1-extra1,1);
stop1  = min(x1+extra1,size(img,1));
start2 = max(x2-extra2,1);
stop2  = min(x2+extra2,size(img,2));
img_region = img(start1:stop1, start2:stop2, :);
spectrum = squeeze(mean(mean(img_region,2),1))';
patch([start2; stop2; stop2; start2],[start1; start1; stop1; stop1],'k')


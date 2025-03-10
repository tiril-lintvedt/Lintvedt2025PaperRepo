function img = find_highfat_pixels(img,threshold, plot_it, mode)
% Finds and selects the high fat pixels in an image structure img, setting 
% rest of image to zero. The selection is based on thresholding the ratio 
% between fat associated peak at 1723 nm and the non peak pixel at 1745 nm.
% mode = 'reflectance' or 'absorbance'
%
%    ratio = relevantPixelsSnv(1723 nm)/ relevantPixelsSnv(1745 nm); 
%    highfat pixels = ratio > threshold
%
% For use on reflectance data (Spectra are "upside down")
    
    
    % Vectorize/flatten the cube (one row correspond to one pixel): 
    %cubeSize = size(img.d);
    mask = mean(img.d,3) ~= 0; % Mask to remove the background (assumed set to 0)
    vectorizedCube = reshape(img.d, [size(img.d,1)*size(img.d,2), size(img.d,3)]);
    vectorizedMask = reshape(mask, [size(img.d,1)*size(img.d,2), 1]);

    relevantPixels = vectorizedCube(vectorizedMask == 1,:); % Non-background pixels
    
    clear mask vectorizedCube

    % SNV on non-background pixels
    %relevantPixelsSnv = snv(relevantPixels);
    
    % Ratio of intensity at I(1723)/I(1745) for each pixel;
    WN1 = 1723;
    WN2 = 1745;
    [~,i1]= min(abs(img.v - WN1)); 
    [~,i2]= min(abs(img.v - WN2)); 
    ratios = relevantPixels(:,i1)./ relevantPixels(:,i2); 
%     mr = mean(ratios);
%     minr = min(ratios);
%     maxr = max(ratios);
    
    % Find high fat pixels and set lowfat pixels to zero
    if strcmp(mode, 'reflectance')
         highfat = ratios < threshold; % "<" because we have reflectance data. Spectra are upside down
    elseif strcmp(mode, 'absorbance')
         highfat = ratios > threshold; 
    end
    relevantPixels(~highfat,:) = zeros(size(relevantPixels(~highfat,:),1), size(img.d,3));
    
    clear ratios
    
%     % Control plot     
%     if plot_it
%         figure
%         plot(img.v,mean(relevantPixels(highfat,:)))
%         title('Mean spectra')
%         text(0,0.95,['Mean-Min-Max ratio: ', str(mr),'-',str(minr),'-',str(maxr)],'Unit','Normalized')
%     end
    
    clear highfat
    
    % Represent the results of the pixel-wise processing in the original
    % spatial arrangement
    vectorizedResults = zeros(size(img.d,1)*size(img.d,2), size(img.d,3)); 
    vectorizedResults(vectorizedMask==1,:) = relevantPixels;
    img.d = reshape(vectorizedResults, [size(img.d,1),size(img.d,2),size(img.d,3)]);
    
   

end
function Xsel = select_channels(X, minBand, maxBand)
% X is a cell array of image structures. This function selects/isloates
% the channels specified by the lower and upper wavelength limits minBand 
% and maxBand in each image. 

   Xsel = {};
   
    for i_img = 1:length(X)

        image = X{i_img};
        cube = image.d;
        Wavelength = image.v;

        [~,minBandIndex] = min(abs(Wavelength-minBand));
        [~,maxBandIndex] = min(abs(Wavelength-maxBand));

        cube = cube(:,:,minBandIndex:maxBandIndex);
        Wavelength = Wavelength(minBandIndex:maxBandIndex);
        
        IMG.d = cube;
        IMG.v = Wavelength;
        IMG.i = image.i;
        Xsel{i_img} = IMG;

    end

end

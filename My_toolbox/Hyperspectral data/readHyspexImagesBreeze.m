function X = readHyspexImagesBreeze(path)
% Read all images in the path to a single cell structure X. Works for 
% images exported in matlab format from Breeze (Prediktor software)
%
% Input:
%           path - string with path to map containing all images
% Output:
%           X - cell array of structures (One structure is one image)

    X = {};
    files = dir(path); files = files(3:end,:);
    [~,idx] = sort([files.datenum]);
    files = files(idx);

    for i_file = 1:length(files)

        % 1) Load the data and the wavelengths (From Breeze export)
        fileName = [path,'\',files(i_file).name,'\',files(i_file).name,'.mat'];
        img = load(fileName); fieldname = fieldnames(img);
        cube = img.(fieldname{1});

        fileName = [path,'\',files(i_file).name,'\',files(i_file).name,'_layers.mat'];
        load(fileName)
        Wavelength = str2num(Wavelength);

        IMG.d = cube;
        IMG.v = Wavelength;
        IMG.i = files(i_file).name;
        X{i_file} = IMG;

    end


end
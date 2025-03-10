function [hyIm,wl] = readHyspexImage(hyspexFile)
% READHYSPEXIMAGE - read ENVI format hyperspectral image from HySpex camera
%
%   Usage:
%   [hyIm,wl] = readHyspexImage(file)
%
%   Input parameters:
%   hyspexFile  -   path to *.hyspex file to be read. Header file *.hdr 
%                   with same name is assumed to be present in same folder.
%               
%   Output parameters:
%   hyIm    -   hyperspectral image, dimensions [nLines, nSamples, nBands]
%   wl      -   wavelengths, vector with nBands elements
%
%   The function operates by extracting image format information from the 
%   header file using regexp, and reading the binary image using 
%   multibandread.
%
%   2016-03-16  Martin H. Skjelvareid, mhs@nofima.no
%
%   Edits:
%       2018-08-28  Stein-Kato Lindberg, stein-kato.lindberg@nofima.no
%       -   Added support for double-precision input image. 
%
%       2018-08-30  Stein-Kato Lindberg, stein-kato.lindberg@nofima.no
%       -   Added support for various data types. 
%       -   Fixed bug which prevented the header file from being read from
%           the current directory. The fileread function returns an error
%           if the input path starts with \. This is not caught by the exist
%           function, which returns 2 in either case (provided the file
%           exists on the path). 

%% Construct header file name and validate paths
[pathstr, name] = fileparts(hyspexFile);

if isempty(pathstr)
    hdrFile = [name '.hdr'];
else
    hdrFile = [pathstr filesep() name '.hdr'];
end

if not(exist(hyspexFile,'file'))
    error(['Hyspex file not found: ' hyspexFile])
end

if not(exist(hdrFile,'file'))
    error(['Header file not found: ' hdrFile])
end

%% Read header file
text = fileread(hdrFile);

%% Find image parameters
tmp = regexp(text,'lines = \d+\s','match');
nLines = str2double(tmp{1}(9:(end-1)));

tmp = regexp(text,'bands = \d+\s','match');
nBands = str2double(tmp{1}(9:(end-1)));

tmp = regexp(text,'samples = \d+\s','match');
nSamples = str2double(tmp{1}(11:(end-1)));

tmp = regexp(text,'header offset = \d+\s','match');
hdrOffset = str2double(tmp{1}(17:(end-1)));

tmp = regexp(text,'interleave = \w+\s','match');
interleave = tmp{1}(14:(end-1));

tmp = regexp(text,'byte order = \d+\s','match');
bigEndian = str2double(tmp{1}(14:(end-1)));        

tmp = regexp(text,'data type = \d+\s', 'match');
dataTypeNumber = str2double(tmp{1}(12:(end-1)));
switch dataTypeNumber
    case 1
        dataType = 'uint8';
    case 2
        dataType = 'int16';
    case 3
        dataType = 'int32';
    case 4
        dataType = 'single';
    case 5
        dataType = 'double';
    case 12
        dataType = 'uint16';
    case 13
        dataType = 'uint32';
    case 14
        dataType = 'int64';
    case 15
        dataType = 'uint64';
    otherwise 
        error('Invalid data type, check the header file. ')
end

%% Find wavelengths
pat1 = 'wavelength = {[\d,\. ]+}';          % Find wavelength text block
tmp = regexp(text,pat1,'match');
wlText = tmp{1}; 
pat2 = '\d+\.\d+';                          % Find individual wavelengths
tmp = regexp(wlText,pat2,'match');
wl = double(length(tmp));
for ii = 1:length(tmp)
    wl(ii) = str2double(tmp{ii});
end

%% Set byte order string for multibandread
if bigEndian                                        
    byteOrder = 'ieee-be';                          
else
    byteOrder = 'ieee-le';
end

%% Read binary file
hyIm = multibandread(hyspexFile,[nLines,nSamples,nBands],dataType,...
    hdrOffset,interleave,byteOrder);

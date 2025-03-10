function [ZQuality]=quality_test(ZSaisir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%  Achim Kohler                                                                        %
%  Center for Biospectroscopy and Data Modelling                                       %
%  Nofima Mat                                                                          %
%  Norwegian Food Research Institute                                                   %
%  Osloveien 1                                                                         %
%  1430 Ås                                                                             %
%  Norway                                                                              %
%                                                                                      %
%  First version: 12.03.05                                                             %
%  Revised: 06.07.21 (TAL)                                                                   %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [QTmatrix]=quality_test(ZSaisir);                                          %
%                                                                                      %
%  Runs a quality test on a Saisir structure                                           %
%  Literature: See Bruker user manual for microorganisms                               % 
%                                                                                      %
%                                                                                      %
%  Status: Running                                                                     %
%                                                                                      %
%  Input:   Saisir structure                                                           %
%           option (0: work with original spectra , 1: first derivative)               %                   
%                                                                                      %
%  Output:  
% Absorbance: max-min absorbance in range 2100-1600
% AmideI : max-min in range 1 (1700-1600), based on first derivative
% Poly: max-min in range 2 (1200-960), based on first derivative
% Noise: max-min in area 2100-2000, based on first derivative
% AbsNoise: std of diff.between raw and smoothed spectrum (Absolute noise)
% AmideIN: AmideI/Noise
% AmideIN2: AmideI/AbsNoise
% PolyN: Poly/Noise
% Water Vapor: max-min in range 1847-1837, based on first derivative
% AmideW: AmideI/Water vapour
% PolyW: Poly/Water vapour
% Ester: max-min in region 1800-1700 (fat content), based on raw spectra
% Fringes: max-min in region 2200-2000, based on raw spectra
% Phos: max-min in region 900-1030, based on raw spectra
% PhosN2: Phos/AbsNoise
% Olefin: max-min in region 1600-1720, based on raw spectra 
% OlefinN2: Olefin/AbsNoise 
% AvgIntN2: Average spectrum intensity/AbsNoise
%

% Recomended cutoffs for removing spectra:
% Absorbance >1.5 or <0.1
% AmideIN<20
% PolyN<10
% AmidIW<20
% PolyW<4
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
ZSaisir1Deriv = saisir_derivative(ZSaisir,2,61,1); % 1st derivative spectrum
ZSaisirSmooth = saisir_derivative(ZSaisir,2,9,0); % Strongly smoothed spectrum 

% Identify the wavenumber range of input data to determine which quality checks
% are possible for the given data:
WNLowerLimit = str2num(ZSaisir.v(1,:));
WNUpperLimit = str2num(ZSaisir.v(end,:));

% Ranges used in subsequent tests:
AbsRange = [1600, 2100]; if WNLowerLimit > min(AbsRange) || WNUpperLimit < max(AbsRange), OutOfRangeAbs = 1; else, OutOfRangeAbs = 0;  end
AmideIRange = [1600, 1700]; if WNLowerLimit > min(AmideIRange) || WNUpperLimit < max(AmideIRange), OutOfRangeAmideI = 1; else, OutOfRangeAmideI = 0; end
PolyRange = [960, 1200]; if WNLowerLimit > min(PolyRange) || WNUpperLimit < max(PolyRange), OutOfRangePoly = 1; else, OutOfRangePoly = 0;end
NoiseRange = [2000, 2100]; if WNLowerLimit > min(NoiseRange) || WNUpperLimit < max(NoiseRange), OutOfRangeNoise = 1; else, OutOfRangeNoise = 0;end
Noise2Range = []; % Can be used for all ranges
WaterRange = [1837, 1847]; if WNLowerLimit > min(WaterRange) || WNUpperLimit < max(WaterRange), OutOfRangeWater = 1; else, OutOfRangeWater = 0;end
EsterRange = [1700, 1800]; if WNLowerLimit > min(EsterRange) || WNUpperLimit < max(EsterRange), OutOfRangeEster = 1;else, OutOfRangeEster = 0; end
FringesRange = [2000, 2200]; if WNLowerLimit > min(FringesRange) || WNUpperLimit < max(FringesRange), OutOfRangeFringes = 1; else, OutOfRangeFringes = 0; end
PhosRange = [900, 1030]; if WNLowerLimit > min(PhosRange) || WNUpperLimit < max(PhosRange), OutOfRangePhos = 1; else, OutOfRangePhos = 0;end
OlefinRange = [1600, 1720]; if WNLowerLimit > min(OlefinRange) || WNUpperLimit < max(OlefinRange), OutOfRangeOlefin = 1; else, OutOfRangeOlefin = 0; end

% -------------------------------------------------------------------------
% QUALITY TEST
% -------------------------------------------------------------------------
% average intensity of each spectrum
AverageIntensity = mean(ZSaisir.d,2);

% absorbance 1600-2100, here OPUS uses raw spectra
WN1=1600.0;
WN2=2100.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir1600_2100=selectcol(ZSaisir,[i1:i2]);
[MaxAbs]=max(ZSaisir1600_2100.d,[],2);
[MinAbs]=min(ZSaisir1600_2100.d,[],2);

% (Signal 1) Amide 1: 1600-1700 (according to OPUS) - First derivative
WN1=1600.0;
WN2=1700.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir1600_1700=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxAmide1]=max(ZSaisir1600_1700.d,[],2);
[MinAmide1]=min(ZSaisir1600_1700.d,[],2);

% (Signal 2) Region 960-1200 - First derivative
WN1=960.0;
WN2=1200.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir960_1200=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxPoly]=max(ZSaisir960_1200.d,[],2);
[MinPoly]=min(ZSaisir960_1200.d,[],2);
 
% Noise: 2000-2100 (according to OPUS) - First derivative
WN1=2000.0;
WN2=2100.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir2000_2100=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxNoise]=max(ZSaisir2000_2100.d,[],2);
[MinNoise]=min(ZSaisir2000_2100.d,[],2);

% Absolute noise (standard deviation - Alternative method based on Guo et al 2020) 
Noise2.d  = ZSaisir.d - ZSaisirSmooth.d; % Isolate noise
Noise2.d = Noise2.d(:,70:(end-70)); % shave off ends to avoid SG edge effects
Noise2.i = ZSaisir.i;
Noise2.v = ZSaisir.v(70:(end-70),:);
Noise2_val = std(Noise2.d,0,2); % Compute standard deviation along row axis

% include plot to control that noise isolation and end-shave was successful:
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
subplot(3,1,1)
plot(str2num(ZSaisir.v),ZSaisir.d)
xlim([str2num(ZSaisir.v(1,:)) str2num(ZSaisir.v(end,:))])
title('Noise estimation control check- Absolute noise','FontSize',16)
ylabel('Original', 'FontSize',14)
subplot(3,1,2)
plot(str2num(ZSaisirSmooth.v),ZSaisirSmooth.d)
xlim([str2num(ZSaisir.v(1,:)) str2num(ZSaisir.v(end,:))])
ylabel('Strongly smoothed', 'FontSize',14)
subplot(3,1,3)
plot(str2num(ZSaisir.v(70:(end-70),:)), Noise2.d)
xlim([str2num(ZSaisir.v(1,:)) str2num(ZSaisir.v(end,:))])
ylabel('Isolated noise','FontSize',14)
xlabel('Raman shift (cm{^{-1}})', 'FontSize',14)
set(gcf, 'Color', [1 1 1])

% Water vapour: 1837-1847 (according to OPUS)- First derivative
WN1=1837.0;
WN2=1847.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir1837_1847=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxWater]=max(ZSaisir1837_1847.d,[],2);
[MinWater]=min(ZSaisir1837_1847.d,[],2);

% Ester: 1700-1800 (gives a measure for the fat content IN RAW SPECTRA)
WN1=1700.0;
WN2=1800.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir1700_1800=selectcol(ZSaisir,[i1:i2]);
[MaxEster]=max(ZSaisir1700_1800.d,[],2);
[MinEster]=min(ZSaisir1700_1800.d,[],2);

% Fringes 2000-2200, here OPUS raw spectra
WN1=2000.0;
WN2=2200.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir2000_2200=selectcol(ZSaisir,[i1:i2]);
[MaxFringes]=max(ZSaisir2000_2200.d,[],2);
[MinFringes]=min(ZSaisir2000_2200.d,[],2);

% Phosphate (e.g ash): Region 900-1030 - Use RAW spectra
WN1=900.0;
WN2=1030.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir900_1030=selectcol(ZSaisir,[i1:i2]);
[MaxPhos]=max(ZSaisir900_1030.d,[],2);
[MinPhos]=min(ZSaisir900_1030.d,[],2);

% Olefinic stretch (unsaturated fatty acids): Region 1600-1720 - Use RAW
% spectra
WN1=1600.0;
WN2=1720.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir1600_1720=selectcol(ZSaisir,[i1:i2]);
[MaxOlefin]=max(ZSaisir1600_1720.d,[],2);
[MinOlefin]=min(ZSaisir1600_1720.d,[],2);


if ~OutOfRangeAbs
    Abs=MaxAbs-MinAbs ; else, Abs = NaN(size(ZSaisir,1),1) ; end
if ~OutOfRangeAmideI 
    AmideI=MaxAmide1-MinAmide1; else, AmideI = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangePoly
    Poly=MaxPoly-MinPoly ; else, Poly = NaN(size(ZSaisir,1),1) ; end
if ~OutOfRangeNoise 
    Noise=MaxNoise-MinNoise; else, Noise = NaN(size(ZSaisir,1),1) ; end 
Noise2 = Noise2_val; % Standard deviation of absolute noise (self-implemented)
AvgIntN2 = AverageIntensity./Noise2;
if ~OutOfRangeWater 
    Water=MaxWater-MinWater; else, Water = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangeEster 
    Ester=MaxEster-MinEster; else, Ester = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangeFringes
    Fringes=MaxFringes-MinFringes; else, Fringes = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangePhos
    Phos = MaxPhos - MinPhos; else, Phos = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangeOlefin
    Olefin = MaxOlefin - MinOlefin; else, Olefin = NaN(size(ZSaisir,1),1) ; end  % Olefinic stretch (unsaturated fatty acids)
if ~OutOfRangeAmideI && ~OutOfRangeNoise 
    AmideIN=AmideI./Noise; else, AmideIN = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangePoly && ~OutOfRangeNoise 
    PolyN=Poly./Noise; else, PolyN = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangeAmideI  
    AmideIN2 =AmideI./Noise2; else, AmideIN2 = NaN(size(ZSaisir,1),1) ; end % Self-implemented
if ~OutOfRangeOlefin  
    OlefinN2 = Olefin./Noise2; else, OlefinN2 = NaN(size(ZSaisir,1),1) ; end % Self-implemented
if ~OutOfRangePhos  
    PhosN2 = Phos./Noise2; else, PhosN2 = NaN(size(ZSaisir,1),1) ; end % Self-implemented
if ~OutOfRangeAmideI && ~OutOfRangeWater 
    AmideW=AmideI./Water; else, AmideW = NaN(size(ZSaisir,1),1) ; end 
if ~OutOfRangePoly && ~OutOfRangeWater 
   PolyW=(MaxPoly-MinPoly)./Water; else, PolyW = NaN(size(ZSaisir,1),1) ; end 

% IdentLowNoise= find(Noise<0.00000000001);
% AmideIN(IdentLowNoise)=0.0;
% PolyN(IdentLowNoise)=0.0;

% IdentLowWater=find(Water<0.00000000001);
% AmideW(IdentLowWater)=0.0;
% PolyW(IdentLowWater)=0.0;

[Nx Ny]=size(ZSaisir.d);
ZQuality.d=zeros(Nx,9);
ZQuality.d(:,1)=Abs;
ZQuality.d(:,2)=AmideI;
ZQuality.d(:,3)=Poly;
ZQuality.d(:,4)=Noise;
ZQuality.d(:,5)=Noise2;
ZQuality.d(:,6)=AmideIN;
ZQuality.d(:,7)=AmideIN2;
ZQuality.d(:,8)=PolyN;
ZQuality.d(:,9)=Water;
ZQuality.d(:,10)=AmideW;
ZQuality.d(:,11)=PolyW;
ZQuality.d(:,12)=Ester;
ZQuality.d(:,13)=Fringes;
ZQuality.d(:,14)= Phos;
ZQuality.d(:,15)= PhosN2;
ZQuality.d(:,16)= Olefin;
ZQuality.d(:,17)= OlefinN2;
ZQuality.d(:,18)= AvgIntN2;

ZQuality.i=ZSaisir.i;
ZQuality.v=['Absorbance '
            'AmideI     '
            'Poly       '
            'Noise      '
            'AbsNoise   '
            'AmideIN    '
            'AmideIN2   '
            'PolyN      '
            'Water Vapor'
            'AmideW     '
            'PolyW      '
            'Ester      '
            'Fringes    '
            'Phos       '
            'PhosN2     '
            'Olefin     '
            'OlefinN2   '
            'AvgIntN2   '];

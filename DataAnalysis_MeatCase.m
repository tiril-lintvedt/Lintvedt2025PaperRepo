% =========================================================================
% PAPER IV - NIR vs Raman robustness comparison
% =========================================================================
%
% ------------------------------------------------------------------------
% 
% Toolboxes required: emsc, GSTools, My_toolbox, Violinplot-Matlab-master
%
%
% Author: Tiril Lintvedt
% Nofima, Raw materials and process optimisation
% email address: tiril.lintvedt@nofima.no
% Website: 
% March 2025; 
% MATLAB version: R2023b
% OS: Windows 11 Home

%------------- BEGIN CODE --------------

%% SETTINGS 
% ----------------------------------------------------------------------

% Main plot colors
plot_colors =[0 0.5 0.9   ; 
              0.9 0.5 0   ; 
              0.2 0.5 0.2 ; 
              0.9 0.3 0.1 ;
              0.7 0.7 0.2 ; 
              0   0   0   ;
              0.5 0.5 0.5 ;
              0.2 0.7 0.7    ];

set(0,'defaultfigurecolor',[1 1 1])

addpath(genpath('Data'))
addpath(genpath('emsc'))
addpath(genpath('GSTools'))
addpath(genpath('My_toolbox'))
addpath(genpath('Violinplot-Matlab-master'))



%% DATA IMPORT
% -------------------------------------------------------------------------

% Main set data (from Petter) ---------------------------------------------
Xkaiser = load_unscrmat_2_saisir('Data\Meat\Kaiser.mat');
Xmm = load_unscrmat_2_saisir('Data\Meat\MarqMetrix.mat');
Xmnir = load_unscrmat_2_saisir('Data\Meat\MicroNIR.mat');

% Generalize spectrum ids to match each other and reference names
Xkaiser.i = replace(cellstr(Xkaiser.i),'_KA','');
Xkaiser.i = replace(cellstr(Xkaiser.i),'_Kai','');
Xkaiser.i = replace(cellstr(Xkaiser.i),'Fin','Fin0');
Xkaiser.i = char(replace(cellstr(Xkaiser.i),'Hom','Hom0'));

Xmm.i = replace(cellstr(Xmm.i),'_MM','');
Xmm.i = replace(cellstr(Xmm.i),'Fin','Fin0');
Xmm.i = char(replace(cellstr(Xmm.i),'Hom','Hom0'));

Xmnir.i = replace(cellstr(Xmnir.i),'-coarse-','_Grov_');
Xmnir.i = replace(cellstr(Xmnir.i),'-fine-','_Fin0_');
Xmnir.i = replace(cellstr(Xmnir.i),'-hom-','_Hom0_');
Xmnir.i = char(replace(cellstr(Xmnir.i),'L','Y')); 
Xmnir.v = char(replace(cellstr(Xmnir.v),',','.'));

X_instruments = {Xkaiser, Xmm, Xmnir};

% New Flatbiff data -------------------------------------------------------
% -------------------------------------------------------------------------

% Kaiser
spec = GSImportspec('Data\Meat\Flatbiff - new experiment\Kaiser\Spectra',0);
[XkaiserF, XmetaKaiserF] = extractSPC_KA(spec); % Spectra and Meta data
log = readcell('Data\Meat\Flatbiff - new experiment\Kaiser\Flatbiff name log.txt');
names = char(log(:,3));

% Rename X and XMeta according to name log file (Assume correct file order)
XkaiserF.i = names;
XmetaKaiserF.i = names;

XkaiserF.i = char(replace(cellstr(XkaiserF.i),'_Kai',''));
XkaiserF.i = char(replace(cellstr(XkaiserF.i),'_a','a'));
XkaiserF.i = char(replace(cellstr(XkaiserF.i),'_b','b'));
XkaiserF.i = char(replace(cellstr(XkaiserF.i),'_c','c'));

% MM 
spec = GSImportspec('Data\Meat\Flatbiff - new experiment\MM',0);
[XmmF, XmetaMMF] = extractSPC(spec); % Spectra and Meta data
% Note : seem to loose some decimals along during the extractSpec
% function..
XmmF.i = replace(cellstr(XmmF.i),'_MM','');
XmmF.i = char(replace(cellstr(XmmF.i),'_00000',''));
XmmF.i = char(replace(cellstr(XmmF.i),'_a','a'));
XmmF.i = char(replace(cellstr(XmmF.i),'_b','b'));
XmmF.i = char(replace(cellstr(XmmF.i),'_c','c'));

% MikroNIR
load('Data\Meat\Flatbiff - new experiment\Flatbiff_MikroNIR.mat')
XmnirF.d = DataMatrix(:,1:125); 
XmnirF.i = ObjLabels; 
XmnirF.v = char(replace(cellstr(VarLabels0),',','.'));
XmnirF.v = XmnirF.v(1:125,:); 
XmnirF.i = replace(cellstr(XmnirF.i),'-Mikro-','_');
XmnirF.i = replace(cellstr(XmnirF.i),'-','_');
XmnirF.i = replace(cellstr(XmnirF.i),'Fin','Fin0');
XmnirF.i = char(replace(cellstr(XmnirF.i),'Hom','Hom0'));
XmnirF.i = char(replace(cellstr(XmnirF.i),'_a','a'));
XmnirF.i = char(replace(cellstr(XmnirF.i),'_b','b'));
XmnirF.i = char(replace(cellstr(XmnirF.i),'_c','c'));

% Reference data -----------------------------------------------------------
Y = readcell('Data\Meat\Referanse_All.xlsx','Sheet', 'Samlet');

sample_names = Y(3:62,1:3);
nsamp = size(sample_names,1);

References.d = cell2mat(Y(3:62,20:22));
References.v = char(Y(2,20:22));
References.i = [char(sample_names(:,1)),repelem('_',nsamp,1), ...
                char(sample_names(:,2)),repelem('_',nsamp,1), ...
                num2str([sample_names{:,3}]')]; 

% Generalize reference id to match spectrum names
References.i = char(replace(cellstr(References.i),' ','0')); 
References.i = char(replaceBetween(cellstr(References.i),1,10,cellstr(References.i(:,1))));


% Flatbif references with "wet" nmr method ------------

Y = readcell('Data\Meat\Fett_vann_F_NIRvsRaman_nov23_newrefF_wet_methods.xlsx','Sheet', 'Sheet1'); % New ref analyses (wet method)

sample_names = Y(34:63, 3);
nsamp = size(sample_names,1);

ReferencesF.d = cell2mat(Y(34:63,19:21));
ReferencesF.v = char(Y(2,19:21));
ReferencesF.i = char(sample_names);

% Generalize reference id to match spectrum names
ReferencesF.i = char(replace(cellstr(ReferencesF.i),'hom','Hom0_')); 
ReferencesF.i = char(replace(cellstr(ReferencesF.i),'fin','Fin0_')); 
ReferencesF.i = char(replace(cellstr(ReferencesF.i),'grov','Grov_')); 

%% EXPLORE RAW SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser.v), Xkaiser.d);
plot_Raman(str2num(Xmm.v), Xmm.d);
plot_NIR(str2num(Xmnir.v), Xmnir.d) ;

plot_Raman(str2num(XkaiserF.v), XkaiserF.d);
plot_Raman(str2num(XmmF.v), XmmF.d);
plot_NIR(str2num(XmnirF.v), XmnirF.d);


plot_from_identifier(XmmF,3:6)

%% PREPROCESSING - PART I
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xkaiser.v)-515));  % 520
[~,i2] = min(abs(str2num(Xkaiser.v)-1870)); % 1800 % 1770
Xkaiser = selectcol(Xkaiser,i1:i2);               

% MM ----------------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xmm.v)-515));  % 520
[~,i2] = min(abs(str2num(Xmm.v)-1870)); % 1800 % 1770
Xmm = selectcol(Xmm,i1:i2);

% MikroNIR ----------------------------------------------------------------
% Keep full region


% NEW FLATBIFF DATA -------------------------------------------------------
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

% Remove non-meat spectra
XkaiserF_cyclo = select_from_identifier(XkaiserF,1,'Cyclo');
XkaiserF = delete_from_identifier(XkaiserF,1,'Cyclo');
XkaiserF = delete_from_identifier(XkaiserF,1,'alufoil');
XkaiserF = delete_from_identifier(XkaiserF,1,'metalplate');
XkaiserF = delete_from_identifier(XkaiserF,13,'discard');
XkaiserF = delete_from_identifier(XkaiserF,1,'F_Fin0_04_b'); % Faulty scan (one extra scan c added here)
XkaiserF.i = char(replace(cellstr(XkaiserF.i),'F_Fin0_04_c','F_Fin0_04_b'));

% Flip spectrum
XkaiserF.d = fliplr(XkaiserF.d);
XkaiserF.v = num2str(fliplr(str2num(XkaiserF.v)')');

% Raman shift range
[~,i1] = min(abs(str2num(XkaiserF.v)-515));  % 520
[~,i2] = min(abs(str2num(XkaiserF.v)-1870)); % 1800
XkaiserF = selectcol(XkaiserF,i1:i2);

% MM ----------------------------------------------------------------------
XmmF_cyclo = select_from_identifier(XmmF,1,'Cyclo');
XmmF = delete_from_identifier(XmmF,1,'Cyclo');

% Raman shift range
[~,i1] = min(abs(str2num(XmmF.v)-515));  % 520
[~,i2] = min(abs(str2num(XmmF.v)-1870)); % 1800
XmmF = selectcol(XmmF,i1:i2);

% MikroNIR ----------------------------------------------------------------



clear i1 i2



%% EXPLORE NEW SPECTRAL REGIONS
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser.v), Xkaiser.d);
plot_Raman(str2num(Xmm.v), Xmm.d);
plot_NIR(str2num(Xmnir.v), Xmnir.d); 

plot_Raman(str2num(XkaiserF.v), XkaiserF.d);
plot_Raman(str2num(XmmF.v), XmmF.d);
plot_NIR(str2num(XmnirF.v), XmnirF.d); 

%% SPIKE REMOVAL FOR RAMAN SPECTRA 
%  ------------------------------------------------------------------------
% Only needed for MM

% Remove damaged pixels/spikes MM
dpix =[660:664, 950:954, 1140:1144, 1584:1588, 1768:1800, 1980:1983, 1996:1999]; 
[~,di] = min(abs(str2num(Xmm.v)-dpix)); % 1800
Xmm = deletecol(Xmm,di);
th = 15; % threshold for spike detection
[Xmm, ~] = spikefix_whitaker_multi(Xmm,2,2,1,th,1); % double check this later

% Remove damaged pixels/spikes Kaiser
dpix = [1768:1792]; 
[~,di] = min(abs(str2num(Xkaiser.v)-dpix)); % 1800
Xkaiser = deletecol(Xkaiser,di);
XkaiserF =deletecol(XkaiserF,di);
th = 45; % threshold for spike detection
[Xkaiser, ~] = spikefix_whitaker_multi(Xkaiser,2,2,1,th,1); % double check this later


% NEW FLATBIFF DATA -------------------------------------------------------
% -------------------------------------------------------------------------

% MM ----------------------------------------------------------------------
dpix =[660:664, 950:954, 1140:1144, 1584:1588, 1768:1800, 1980:1983, 1996:1999]; 
[~,di] = min(abs(str2num(XmmF.v)-dpix)); % 1800
XmmF = deletecol(XmmF,di);


%% EXPLORE SPIKE FREE SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser.v), Xkaiser.d);
plot_Raman(str2num(Xmm.v), Xmm.d);
plot_NIR(str2num(Xmnir.v), Xmnir.d) ;

plot_Raman(str2num(XkaiserF.v), XkaiserF.d);
plot_Raman(str2num(XmmF.v), XmmF.d);
plot_NIR(str2num(XmnirF.v), XmnirF.d); 


%% SMOOTHING RAMAN SPECTRA
% -------------------------------------------------------------------------
wd = 9;
derorder = 2;

% MM ----------------------------------------------------------------------
Xmm = saisir_derivative(Xmm,derorder,wd,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(Xmm.v);
Xmm = selectcol(Xmm,5:(nvar-4));% Remove edge effects (the four edge points on each side)

% Kaiser ------------------------------------------------------------------
Xkaiser= saisir_derivative(Xkaiser,derorder,wd,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(Xkaiser.v);
Xkaiser = selectcol(Xkaiser,5:(nvar-4));% Remove edge effects (the four edge points on each side)


% NEW FLATBIFF DATA -------------------------------------------------------
% -------------------------------------------------------------------------

% MM ----------------------------------------------------------------------
XmmF = saisir_derivative(XmmF,derorder,wd,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(XmmF.v);
XmmF = selectcol(XmmF,5:(nvar-4));% Remove edge effects (the four edge points on each side)

% Kaiser ------------------------------------------------------------------
XkaiserF= saisir_derivative(XkaiserF,derorder,wd,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(XkaiserF.v);
XkaiserF = selectcol(XkaiserF,5:(nvar-4));% Remove edge effects (the four edge points on each side)



%% EXPLORE SMOOTHED SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser.v), Xkaiser.d);
plot_Raman(str2num(Xmm.v), Xmm.d);
plot_NIR(str2num(Xmnir.v), Xmnir.d); 

plot_Raman(str2num(XkaiserF.v), XkaiserF.d);
plot_Raman(str2num(XmmF.v), XmmF.d);
plot_NIR(str2num(XmnirF.v), XmnirF.d); 


%% MERGE OLD AND NEW DATA SETS
% ----------------------------------------------------------------------
% Delete extra replicates
XkaiserF = delete_from_identifier(XkaiserF,10,'c');

% Merge old/new spectral data ----------------
Xkaiser = saisir_rowmerge(Xkaiser, XkaiserF);
Xmm = saisir_rowmerge(Xmm, XmmF);  
Xmnir = saisir_rowmerge(Xmnir, XmnirF); 


% Merge old/new reference data MISSING ------------------------------------
References = saisir_rowmerge(References, ReferencesF);

%%  EXPLORE REFERENCE VALUES
% -------------------------------------------------------------------------

saisir_explore_reference_data(References,3:6, 8:9) % 1


%% SEPARATE INSTRUMENTS AND SUBSETS
% -------------------------------------------------------------------------

% Different region versions
X_instrument_subsets = {Xkaiser, Xmm, Xmnir};
X_subsets = {};
n_subsets = {};

for i = 1:size(X_instrument_subsets,2)
    
    X = X_instrument_subsets{i};
    % Separate by muscle/part
    [X_Ytrefilet, ~] = select_from_identifier(X,1,'Y'); 
    [X_Bankekjott, ~] = select_from_identifier(X,1,'B'); 
    [X_Flatbiff, ~] = select_from_identifier(X,1,'F'); 

    % Separate degree of heterogeneity
    [X_grov, ~] = select_from_identifier(X,3,'Grov');
    [X_fin, ~] = select_from_identifier(X,3,'Fin');
    [X_hom, ~] = select_from_identifier(X,3,'Hom');

    % Gather all data subsets ------------------------------------------------- 
    nsets = 7;
    X_subsets(end+1:end+nsets) = {X_Ytrefilet, X_Bankekjott, X_Flatbiff, X_grov, X_fin, X_hom, X};
    n_subsets(end+1 : end+nsets , 1) = {'Ytrefilet'; 'Bankekjøtt'; 'Flatbiff'; ...
                                    'Grovmalt (20 mm)'; 'Finmalt (4 mm)'; ...
                                    'Homogenisert (Stefan mikser)'; 'ALLE' };

end



%% CONNECT REFERENCE DATA WITH RAMAN SPECTRA
% -------------------------------------------------------------------------
ref_idpos = 1:9; % string position in reference id tags
x_idpos = 1:9; %string position in spectral/data blocks 
reference_names = References.i(:,ref_idpos);
nref = size(References.d,2); % Number of reference analyses
Y_subsets = {};

for i_set = 1:length(X_subsets)  % For each subset, create an R block
    
    Xi = X_subsets{i_set};
    R.d = zeros(size(Xi.i,1),nref);
    R.v = References.v;
    R.i = {};

    for i_sample = 1:size(Xi.i,1)
        name = Xi.i(i_sample,x_idpos);

        for ref = 1:size(reference_names,1)
            ref_name = reference_names(ref,:);       
            if strcmp(name,ref_name)
                R.d(i_sample,:) = References.d(ref,:);
                R.i{i_sample} = Xi.i(i_sample,:);
                break
            end     
        end

    end
    R.i = char(R.i');
    Y_subsets{i_set} = R;

end

 
% Save for Unscrambler
Ykaiser_all = Y_subsets{7};  
Ymm_all = Y_subsets{14}; 
Ymnir_all = Y_subsets{21};  

save('Data\ReferencesAll.mat','Ykaiser_all','Ymm_all','Ymnir_all')  % Matlab saisir version

clear Ykaiser_all Ymm_all Ymnir_all x_idpos Xi i_set i_sample name nref ...
      ref_idpos ref ref_name refernce_names


%% BASELINE CORRECTION ONLY (ALS) FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

X_prep_als_subsets = {};
Fbaseline_subsets = {};
lambda = 4.5; % Smoothing parameter
p = 0.001;

for i = 1:7 % Kaiser
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    X_prep_als_subsets{i} = Subset_prep;

    Baselines.i = Subset.i;
    Baselines.v = Subset.v;
    Baselines.d = baseline;
    Fbaseline_subsets{i} = Baselines;
    
    % Control plot
    plot_Raman(str2num(Subset_prep.v), Subset_prep.d);

end

% MM ----------------------------------------------------------------------

lambda = 3.5; % Smoothing parameter 
p = 0.001;
for i = 8:14 % MM
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    X_prep_als_subsets{i} = Subset_prep;
    
    Baselines.i = Subset.i;
    Baselines.v = Subset.v;
    Baselines.d = baseline;
    Fbaseline_subsets{i} = Baselines;

    % Control plot
    plot_Raman(str2num(Subset_prep.v), Subset_prep.d);

end


%% SNV FOR NIR SPECTRA
% -------------------------------------------------------------------------

X_prep_snv_subsets = {};

for i = 15:21
    Subset = X_subsets{i};
    [Subset_prep] = saisir_snv(Subset); % SNV correction
    X_prep_snv_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);

end


%% EXPLORE PREP NIR SPECTRA
% -------------------------------------------------------------------------
isubset = 21;
plot_NIR(str2num(X_prep_snv_subsets{isubset}.v), X_prep_snv_subsets{isubset}.d); ylim([-2 2]); grid on

fatIndicator = str2num(X_prep_snv_subsets{isubset}.i(:,8:9));
fatPercent = Y_subsets{isubset}.d(:,1);
color_by(fatPercent)
cb = colorbar;
title(cb,'%Fat', 'Fontsize', 14)

%% SNV FOR NIR SPECTRA, FAT BAND ALONE
% -------------------------------------------------------------------------

X_prep_snv_fatband_subsets = {};

for i = 15:21
    Subset = X_subsets{i};
    Subset = selectwn(Subset,1149:1280); % SG 
    [Subset_prep] = saisir_snv(Subset); % SNV correction
    X_prep_snv_fatband_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);

end



%% FAT MODEL FOR TEST ON OXYGENATION SET
% -------------------------------------------------------------------------

% MM ------------------------------------------------------------------

X = X_prep_als_subsets{14}; % Spectra
Y = Y_subsets{14}; % Reference values


% Look only at one grinding type only
type = 'Hom'; % 'B','F','Y'
X_type = select_from_identifier(X,3:5,type); 
Y_type = select_from_identifier(Y,3:5,type);

ncomp = 1;
target = 1;
[b0,B,T,W,P,q] = pls_nipals(X_type.d,Y_type.d(:,target),ncomp);
YpredCal.d = b0 + X_type.d*B;

%plot_regcoef(X.v, B,  0.5,'findPeaks', 'no');
plot(Y_type.d(:,1),YpredCal.d,'O')

TotFatModel.B = B;
TotFatModel.b0 = b0;
TotalFatModel.v = X_type.v;
TotFatModel.i =  'Prediction by: Ypred = b0 + X*B';


%% OVERALL VALIDATION PERFORMANCES - LOOCV GENERAL MODELS
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------
X = X_prep_als_subsets{7}; % Spectra
X = selectwn(X,524:1770); % Remove tail with broad interferant peak
Y = Y_subsets{7}; % Reference values

% Model frames
ncomp = 15; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 3:9; % idpos defining replicates to be held out simultaneously in CV for aopt choice.

[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X, Y, target, reps_idpos, ncomp, 'plotLimits', [-3 25]);


% Plot the regression coefficients
plot_regcoef(Regcoeff.v, squeeze(Regcoeff.d(:,:,Performance.lv)),0.5,'scaledPlot','no','findPeaks','no')
legend(cellstr(Regcoeff.i),'Location','northwest')
xlabel('Raman shift (cm^{-1})', 'FontSize',18)

% MM ------------------------------------------------------------------

X = X_prep_als_subsets{14}; % Spectra
X = selectwn(X,524:1770); % Remove tail with broad interferant peak
Y = Y_subsets{14}; % Reference values


% Model frames
ncomp = 15; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 3:9; % idpos defining replicates to be held out simultaneously in CV for aopt choice.

[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X, Y, target, reps_idpos, ncomp, 'plotLimits', [-3 25]);


% Plot the regression coefficients
plot_regcoef(Regcoeff.v, squeeze(Regcoeff.d(:,:,Performance.lv)),0.5,'scaledPlot','no','findPeaks','no')
legend(cellstr(Regcoeff.i),'Location','northwest')
xlabel('Raman shift (cm^{-1})', 'FontSize',18)


% MikroNIR ------------------------------------------------------------------

X = X_prep_snv_fatband_subsets{21}; % Spectra
Y = Y_subsets{21}; % Reference values


% Model frames
ncomp = 15; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 3:9; % idpos defining replicates to be held out simultaneously in CV for aopt choice.

[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X, Y, target, reps_idpos, ncomp, 'plotLimits', [-3 25]);

% Plot the regression coefficients
plot_regcoef(Regcoeff.v, squeeze(Regcoeff.d(:,:,Performance.lv)),0.5,'scaledPlot','no','findPeaks','no')
legend(cellstr(Regcoeff.i),'Location','northwest')
xlabel('Wavelength (nm)', 'FontSize',18)


%% ALL COMBINATIONS VALIDATION SCHEME 
% -------------------------------------------------------------------------

% KAISER ------------------------------------------------------------------
X = X_prep_als_subsets{7}; % Spectra
X = selectwn(X,524:1770); % Remove tail with broad interferant peak
Y = Y_subsets{7}; % Reference values

% Model frames
aopt = 1; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 1:9; % idpos defining replicates to be held out simultaneously in CV for aopt choice.
[FIG, Regcoeff, Performance, Residuals] = PLSR_allcombos_validation_from_identifier(X, Y, target, idpos, aopt, reps_idpos, 'aoptAlg','STATIC', 'aoptAlgMetric',...
                                            'RMSEP_corr','aoptMode','calCV' , 'plotLimits', [-5 50]);
% Update pred vs target plot 
leg = get(gca, 'Legend');
set(leg,'Location', 'southeast')
newlegstr = strrep(leg.String,'Fin0','Medium');
newlegstr = strrep(newlegstr,'Grov','Coarse');
newlegstr = strrep(newlegstr,'Hom0','Fine');
newlegstr = strrep(newlegstr,'( LV = 1 )','');
newlegstr = strrep(newlegstr,'->',' -> ');
leg.String = newlegstr;
text(0.1,0.95, ['LV = ',num2str(aopt)],'Unit', 'Normalized','Fontsize', 14)
xlim([-5 30])
ylim([-5 30])

% Plot the regression coefficients
plot_regcoef(Regcoeff.v, Regcoeff.d,0.5,'scaledPlot','no','findPeaks','no','plotColors', plot_colors)
xline(1655, 'Color', [0.5 0.85 0.5],'Linewidth', 5) % Indicate fat related peaks (Beattie et al):
xline(1445, 'Color', [0.5 0.85 0.5],'Linewidth', 35)
xline(1301, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1267, 'Color', [0.5 0.85 0.5],'Linewidth', 10)
xline(1085, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1744, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1675, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1128, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1063, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(860, 'Color', [0.5 0.85 0.5],'Linewidth', 40)
xlim([500 1800])
legend(cellstr(Regcoeff.i),'Location','Eastoutside')
xlabel('Raman shift (cm^{-1})', 'FontSize',18)


% Plot performance summary
Performance.i = char(strrep(cellstr(Performance.i),'Hom0','Fine'));
Performance.i = char(strrep(cellstr(Performance.i),'Fin0','Medium'));
Performance.i= char(strrep(cellstr(Performance.i),'Grov','Coarse'));
Performance.i= char(strrep(cellstr(Performance.i),'(Calibration)',''));
Performance.v = char(strrep(cellstr(Performance.v),'Hom0','Fine'));
Performance.v = char(strrep(cellstr(Performance.v),'Fin0','Medium'));
Performance.v= char(strrep(cellstr(Performance.v),'Grov','Coarse'));
Performance.v= char(strrep(cellstr(Performance.v),'(Validation)',''));
[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 4])
set(ah{2}, 'YLim',[-5.5 5.5])
set(ah{3}, 'YLim',[0 1.5])


% Check error correlations between grinding degrees combos  ---------------

% Indicator for Raman scattering intensity  (General Raman scattering level) 
X_bslcorintavg = X_prep_als_subsets{7};
X_bslcorintavg.d = mean(X_bslcorintavg.d,2);

% Indicator for reflection level in dataset (General raw absorbance level NIR)
X_rawabsavg = X_subsets{21};
X_rawabsavg.d = mean(X_rawabsavg.d,2);

% Indicator for fluorescence level in data set (General baseline level)
X_fluorescence = Fbaseline_subsets{7};
X_fluorescence.d = mean(X_fluorescence.d,2);

X_Ramsc_avg = average_from_identifier(X_bslcorintavg,3:6); % Raman scattering intensity collected (Several reasons why intensity could be lower or higher)
X_Ramsc_avg.v = 'Avg Raman scattering';
X_refl_avg = average_from_identifier(X_rawabsavg,3:6); % NIR scattering indication
X_refl_avg.v = 'Avg Raw NIR abs';
X_fl_avg = average_from_identifier(X_fluorescence,3:6);
X_fl_avg.v = 'Avg Fluorescence';
X_fl_std = std_from_identifier(X_fluorescence,3:6);
X_fl_std.v = 'Std Fluorescence';

X_fl_avg.i = ['Medium';'Coarse';'Fine  '];
xcats = categorical(cellstr(X_fl_avg.i));
xcats = reordercats(xcats,{'Fine' 'Medium' 'Coarse'});
figure; hold on
bar(xcats,X_fl_avg.d)
ylabel('Fluorescence indicator', 'Fontsize', 18); box off
grid on
ylabel({'Reflection property indicator';'(Avg. intensity, NIR raw abs)'}, 'Fontsize', 14); box off


% Use differences in scattering indicators between group combos and
% evaluate correlations 
names = {'fin0', 'grov', 'hom'};
IndNames = {'Avg Refl.', 'Avg Ram. intensity','Avg Fluor.','Std Fluor.'};
Performancev2 = Performance;
Performancev2.slope = 1 - Performancev2.slope; 
[DiffScattind_vect, Diffperformance_vect, description] =  check_error_correlations(Performancev2, X_fl_avg, names, 'slope',0, 'showPlot','yes');
plot_error_correlations(Performancev2,description,IndNames, X_refl_avg, X_Ramsc_avg, X_fl_avg, X_fl_std);



% MM ----------------------------------------------------------------------

X = X_prep_als_subsets{14}; % Spectra
X = selectwn(X,524:1770); % Remove tail for correspondence with Kaiser set.
Y = Y_subsets{14}; % Reference values


% Model frames
aopt = 1; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 1:9; % idpos defining replicates to be held out simultaneously in CV for aopt choice.
[FIG, Regcoeff, Performance, Residuals] = PLSR_allcombos_validation_from_identifier(X, Y, target, idpos, aopt, reps_idpos, 'aoptAlg','STATIC', 'aoptAlgMetric',...
                                            'RMSEP_corr','aoptMode','calCV' , 'plotLimits', [-5 50]);
% Update pred vs target plot 
leg = get(gca, 'Legend');
set(leg,'Location', 'southeast')
newlegstr = strrep(leg.String,'Fin0','Medium');
newlegstr = strrep(newlegstr,'Grov','Coarse');
newlegstr = strrep(newlegstr,'Hom0','Fine');
newlegstr = strrep(newlegstr,'( LV = 1 )','');
newlegstr = strrep(newlegstr,'->',' -> ');
leg.String = newlegstr;
text(0.1,0.95, ['LV = ',num2str(aopt)],'Unit', 'Normalized','Fontsize', 14)
xlim([-5 30])
ylim([-5 30])

% Plot the regression coefficients
plot_regcoef(Regcoeff.v, Regcoeff.d,0.5,'scaledPlot','no','findPeaks','no','plotColors', plot_colors)
xline(1655, 'Color', [0.5 0.85 0.5],'Linewidth', 5) % Indicate fat related peaks (Beattie et al):
xline(1445, 'Color', [0.5 0.85 0.5],'Linewidth', 35)
xline(1301, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1267, 'Color', [0.5 0.85 0.5],'Linewidth', 10)
xline(1085, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1744, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1675, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1128, 'Color', [0.5 0.85 0.5],'Linewidth', 5)
xline(1063, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(860, 'Color', [0.5 0.85 0.5],'Linewidth', 40)
legend(cellstr(Regcoeff.i),'Location','Eastoutside')
xlim([500 1800])
xlabel('Raman shift (cm^{-1})', 'FontSize',18)

% Plot performance summary
Performance.i = char(strrep(cellstr(Performance.i),'Hom0','Fine'));
Performance.i = char(strrep(cellstr(Performance.i),'Fin0','Medium'));
Performance.i= char(strrep(cellstr(Performance.i),'Grov','Coarse'));
Performance.i= char(strrep(cellstr(Performance.i),'(Calibration)',''));
Performance.v = char(strrep(cellstr(Performance.v),'Hom0','Fine'));
Performance.v = char(strrep(cellstr(Performance.v),'Fin0','Medium'));
Performance.v= char(strrep(cellstr(Performance.v),'Grov','Coarse'));
Performance.v= char(strrep(cellstr(Performance.v),'(Validation)',''));

[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 4])
set(ah{2}, 'YLim',[-5.5 5.5])
set(ah{3}, 'YLim',[0 1.5])

% Check correlation between spectrum intensities and reg coef intensity
RegcoefsAvg = average_from_identifier(Regcoeff,1:4);
GrindingDegAvg = average_from_identifier(X,3:6);
corRegInt = corrcoef(mean(RegcoefsAvg.d,2), mean(GrindingDegAvg.d,2));

% Check error correlations across grinding degree combos ------------------

% Indicator for Raman scattering intensity (General Raman scattering level)
X_bslcorintavg = X_prep_als_subsets{14};
X_bslcorintavg.d = mean(X_bslcorintavg.d,2);

% Indicator 2 for scattering level in dataset (General raw absorbance level NIR)
X_rawabsavg = X_subsets{21};
X_rawabsavg.d = mean(X_rawabsavg.d,2);

% Indicator for fluorescence level in data set (General baseline level)
X_fluorescence = Fbaseline_subsets{14};
X_fluorescence.d = mean(X_fluorescence.d,2);

X_Ramsc_avg = average_from_identifier(X_bslcorintavg,3:6); % Raman scattering intensity collected (Several reasons why intensity could be lower or higher)
X_refl_avg = average_from_identifier(X_rawabsavg,3:6); % NIR scattering indication
X_fl_avg = average_from_identifier(X_fluorescence,3:6);
X_fl_std = std_from_identifier(X_fluorescence,3:6);


% Fluorescence indicator plot
X_fl_avg.i = ['Medium';'Coarse';'Fine  '];
xcats = categorical(cellstr(X_fl_avg.i));
xcats = reordercats(xcats,{'Fine' 'Medium' 'Coarse'});
figure; hold on
bar(xcats,X_fl_avg.d)
ylabel('Fluorescence indicator', 'Fontsize', 18); box off
grid on

% Use differences in scattering indicators between group combos
names = ['fin0', 'grov', 'hom'];
Performancev2 = Performance;
Performancev2.slope = 1 - Performancev2.slope; 
[DiffScattind_vect, Diffperformance_vect, description] =  check_error_correlations(Performancev2, X_fl_avg, names, 'slope',0, 'showPlot','no');

% Use differences in scattering indicators between group combos and
% evaluate correlations 
names = {'fin0', 'grov', 'hom'};
IndNames = {'Avg Refl.', 'Avg Ram. intensity','Avg Fluor.','Std Fluor.'};
Performancev2 = Performance;
Performancev2.slope = 1 - Performancev2.slope; 
plot_error_correlations(Performancev2,description,IndNames, X_refl_avg, X_Ramsc_avg, X_fl_avg, X_fl_std);


% MicroNIR --------------------------------------------------------------------


X = X_prep_snv_fatband_subsets{21}; % Spectra
Y = Y_subsets{21}; % Reference values

% Model frames
aopt = 2; % Number of components to use is predefined  
target = 1; % (col position in Y) // Fat only
idpos = 3:6;
reps_idpos = 1:9;
[FIG, Regcoeff, Performance, Residuals] = PLSR_allcombos_validation_from_identifier(X, Y, target, ...
                                            idpos, aopt,reps_idpos, 'aoptAlg','STATIC', 'aoptAlgMetric',...
                                            'RMSEP_corr','aoptMode','calCV' , 'plotLimits', [-5 50]); 
% Update pred vs target plot 
leg = get(gca, 'Legend');
set(leg,'Location', 'southeast')
newlegstr = strrep(leg.String,'Fin0','Medium');
newlegstr = strrep(newlegstr,'Grov','Coarse');
newlegstr = strrep(newlegstr,'Hom0','Fine');
newlegstr = strrep(newlegstr,'( LV = 2 )','');
newlegstr = strrep(newlegstr,'->',' -> ');
leg.String = newlegstr;
text(0.1,0.95, ['LV = ',num2str(aopt)],'Unit', 'Normalized','Fontsize', 14)
xlim([-5 30])
ylim([-5 30])

% Plot the regression coefficients
Regcoeff.i = char(strrep(cellstr(Regcoeff.i),'Hom0','Fine'));
Regcoeff.i = char(strrep(cellstr(Regcoeff.i),'Fin0','Medium'));
Regcoeff.i= char(strrep(cellstr(Regcoeff.i),'Grov','Coarse'));
FIG2 = plot_regcoef(Regcoeff.v, Regcoeff.d,0.5,'scaledPlot','no','findPeaks','no','plotColors', plot_colors);
xline(1211, 'Color', [0.5 0.85 0.5],'Linewidth', 20)
legend(cellstr(Regcoeff.i),'Location', 'Eastoutside')
xlabel('Wavelength', 'FontSize',18)


% Plot performance summary
Performance.i = char(strrep(cellstr(Performance.i),'Hom0','Fine'));
Performance.i = char(strrep(cellstr(Performance.i),'Fin0','Medium'));
Performance.i= char(strrep(cellstr(Performance.i),'Grov','Coarse'));
Performance.i= char(strrep(cellstr(Performance.i),'(Calibration)',''));
Performance.v = char(strrep(cellstr(Performance.v),'Hom0','Fine'));
Performance.v = char(strrep(cellstr(Performance.v),'Fin0','Medium'));
Performance.v= char(strrep(cellstr(Performance.v),'Grov','Coarse'));
Performance.v= char(strrep(cellstr(Performance.v),'(Validation)',''));

[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 4])
set(ah{2}, 'YLim',[-5.5 5.5])
set(ah{3}, 'YLim',[0 1.5])


% Check error correlations across grinding degree combos ------------------

% Indicator for scattering level in dataset (General raw absorbance level)
X_rawabsavg = X_subsets{21}; %X_rawabsavg = selectwn(X_rawabsavg,1149:1280); 
X_rawabsavg.d = mean(X_rawabsavg.d,2);

% Indicator for fluorescence level in data set (General baseline level)
X_fluorescence = Fbaseline_subsets{14}; % 7, 14
X_fluorescence.d = mean(X_fluorescence.d,2);

% Water content in samples
Y_water = selectcol(Y,3);

X_refl_avg = average_from_identifier(X_rawabsavg,3:6); % NIR scattering indication
X_refl_std = std_from_identifier(X_rawabsavg,3:6); % NIR scattering indication
X_fl_avg = average_from_identifier(X_fluorescence,3:6);
X_water_avg = average_from_identifier(Y_water,3:6);


X_refl_avg.i = ['Medium';'Coarse';'Fine  '];
X_refl_avg.d = 1-X_refl_avg.d; % The reflection indicator is actually absorption, thus we take 1-abs = refl
xcats = categorical(cellstr(X_refl_avg.i));
xcats = reordercats(xcats,{'Fine' 'Medium' 'Coarse'});
figure; hold on
bar(xcats,X_refl_avg.d)
ylabel('Reflection indicator', 'Fontsize', 18); box off
grid on
ylim([0 0.2])

% Use differences in scattering indicators between group combos and
% evaluate correlations 
Performancev2 = Performance;
Performancev2.slope = 1 - Performancev2.slope; 
[DiffScattind_vect, Diffperformance_vect, description] =  check_error_correlations(Performancev2, X_fl_avg, names, 'slope',0, 'showPlot','no');
names = {'fin0', 'grov', 'hom'};
IndNames = {'Avg Reflection ind', 'Std Reflection ind'}; 
plot_error_correlations(Performancev2,description,IndNames, X_refl_avg, X_refl_std);
 


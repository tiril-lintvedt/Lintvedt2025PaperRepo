% =========================================================================
% PAPER // Raman and NIR spectroscopy: A discussion of calibration robustness 
% for food quality measurements through two case studies
%
% Case: Strawberries
% =========================================================================
%
% ------------------------------------------------------------------------
% 
% Toolboxes required: Saisir, GSTools, My_toolbox, Violinplot-Matlab-master
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
% -------------------------------------------------------------------------

% Main plot colors

plot_colors =[0 0.5 0.9   ; 
              0.9 0.5 0   ; 
              0.2 0.5 0.2 ; 
              0.9 0.3 0.1 ;
              0.7 0.7 0.2 ; 
              0   0   0   ;
              0.5 0.5 0.5 ;
              0.2 0.7 0.7 ;
              0.1 0.3 0.9 ;
              0.3 0.9 0.3 ;
              0.9 0.9 0.9 ;
              0.1 0.1 0.5 ;];

set(0,'defaultfigurecolor',[1 1 1])

addpath(genpath('Data'))
%addpath(genpath('emsc'))
addpath(genpath('Saisir'))
addpath(genpath('GSTools'))
addpath(genpath('My_toolbox'))
addpath(genpath('Violinplot-Matlab-master'))

%% DATA IMPORT
% -------------------------------------------------------------------------

% Main set data -----------------------------------------------------------
[Xkaiser2021, Xka21month] = load_unscrmat_2_saisir('Data\Strawberries\Kaiser 2021.mat');
[Xkaiser2022, Xka22month] = load_unscrmat_2_saisir('Data\Strawberries\Kaiser 2022.mat');
[Xmnir2021, Xmn21month] = load_unscrmat_2_saisir('Data\Strawberries\MicroNIR 2021.mat');
[Xmnir2022, Xmn22month] = load_unscrmat_2_saisir('Data\Strawberries\MicroNIR 2022.mat');

% Generalize spectrum ids to match each other and reference names
Xkaiser2021.i = [Xkaiser2021.i, repelem('_',size(Xkaiser2021.d,1),1), num2str(Xka21month.d)];
Xkaiser2021.i = char(replace(cellstr(Xkaiser2021.i),'S','')); 
Xkaiser2021.i = char(replace(cellstr(Xkaiser2021.i),' ','')); 
% Make sure all id numbers have 3 digits 
cellnames = cellstr(Xkaiser2021.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
Xkaiser2021.i = char(newcellnames);


Xkaiser2022.i = [Xkaiser2022.i, repelem('_',size(Xkaiser2022.d,1),1), num2str(Xka22month.d)];
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),' ',''));
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'S','')); 
% Make sure all id numbers have 3 digits 
cellnames = cellstr(Xkaiser2022.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d_)',['00',x(1:2)]),cellnames, 'UniformOutput',false);
Xkaiser2022.i = char(newcellnames);
cellnames = cellstr(Xkaiser2022.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
Xkaiser2022.i = char(newcellnames);

Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'MB','')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),' ','')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'-','_')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'.sam','')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'jord_sept21_','')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'.','')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'_A','A')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'A','_A')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'_B','B')); 
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'B','_B')); 
Xmnir2021.v = char(replace(cellstr(Xmnir2021.v),',','.'));
cellnames = cellstr(Xmnir2021.i);
newcellnames = cellfun(@(x) regexprep(x,'(_\d)$',''),cellnames, 'UniformOutput',false);
Xmnir2021.i = char(newcellnames);
Xmnir2021.i = [Xmnir2021.i, repelem('_',size(Xmnir2021.d,1),1), num2str(Xmn21month.d)]; % add month data
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'b','_B'));
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'a','A'));
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),' _','_')); 
Xmnir2021.i(1171,:) = '096_A_1_202109';
Xmnir2021.i(697,:) =  '017_A_1_202109';


Xmnir2022.i = [Xmnir2022.i, repelem('_',size(Xmnir2022.d,1),1), num2str(Xmn22month.d)];
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'jord','')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),' ','')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'-','_')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'murano_','')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'.sam','')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'A','_A'));
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'B','_B')); 
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'a','_A')); 
Xmnir2022.v = char(replace(cellstr(Xmnir2022.v),',','.'));



% Reference data -----------------------------------------------------------

Y = readcell('Data\Strawberries\Ref.xlsx', 'Sheet', 'Sheet1');
sample_names = Y(2:401,[1,13]); 
nsamp = size(sample_names,1);

References.d = [Y{2:401,10}]'; % Brix only
References.v = char(Y(1,10));
References.i = [char(sample_names(:,1)),repelem('_',nsamp,1), ...
                num2str([sample_names{:,2}]')];  

% Generalize reference id to match spectrum names
References.i = char(replace(cellstr(References.i),' ','')); 
References.i = char(replace(cellstr(References.i),'(2)',''));
References.i = char(replace(cellstr(References.i),'S',''));
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d_)',['00',x(1:2)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d\d_\d\d\d\d\d\d)$',[x(1:4),'000_',x(5:10)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);


% Remove spectra without refernce values from spectral blocks -------------
remove_names = {'041_A_1_202109';
                '041_A_2_202109';
                '041_B_1_202109';
                '041_B_2_202109';
                '041_A_3_202109';
                '041_B_3_202109';
                '058_A_1_202109';
                '058_A_2_202109';
                '058_B_1_202109';
                '058_B_2_202109';
                '058_A_3_202109';
                '058_B_3_202109';
                '065_A_1_202109';
                '065_A_2_202109';
                '065_B_1_202109';
                '065_B_2_202109';
                '065_A_3_202109';
                '065_B_3_202109';
                '060_A_1_202106';
                '060_A_2_202106';
                '060_A_3_202106';
                '060_B_1_202106';
                '060_B_2_202106';
                '060_B_3_202106';
                '066_A_1_202106';
                '066_A_2_202106';
                '066_A_3_202106';
                '066_B_1_202106';
                '066_B_2_202106';
                '066_B_3_202106';
                '043_B_2_202106'; % NIR spectra with nan values (replicate still available)
                '068_A_2_202106' }; % NIR spectra with nan values (replicate still available)

for i = 1:length(remove_names)
    Xkaiser2021 = delete_from_identifier(Xkaiser2021,1,remove_names{i,:});
    Xmnir2021 = delete_from_identifier(Xmnir2021,1,remove_names{i,:});
end

% Rename replicates for consistency across sets
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'_3_','_1_'));
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'_4_','_2_'));


%% REMOVE SAMPLES THAT ARE NOT COMMON FOR RAMAN AND NIR SETS
% -------------------------------------------------------------------------

% Remove surplus replicates
Xmnir2021 = delete_from_identifier(Xmnir2021,7,'3');
Xmnir2022 = delete_from_identifier(Xmnir2022,7,'3');


X_subsets_kaiser =  {Xkaiser2021, Xkaiser2022};
X_subsets_mnir = {Xmnir2021, Xmnir2022};

X_subsets_kaiser_intersect = {};
X_subsets_mnir_intersect = {};

for i = 1:length(X_subsets_kaiser)
    
    X_ka_i = X_subsets_kaiser{i};
    X_nir_i = X_subsets_mnir{i};
    
    keep_in_kaiser = find(ismember(cellstr(X_ka_i.i(:,[1:4, 9:14])),cellstr(X_nir_i.i(:,[1:4, 9:14]))));
    keep_in_mnir = find(ismember(cellstr(X_nir_i.i(:,[1:4, 9:14])), cellstr(X_ka_i.i(:,[1:4, 9:14]))));
    
    X_subsets_kaiser_intersect{i} = selectrow(X_ka_i,keep_in_kaiser); 
    X_subsets_mnir_intersect{i} = selectrow(X_nir_i,keep_in_mnir);

    
end

Xkaiser2021 =  X_subsets_kaiser_intersect{1};
Xkaiser2022 = X_subsets_kaiser_intersect{2};
Xmnir2021 = X_subsets_mnir_intersect{1};
Xmnir2022 = X_subsets_mnir_intersect{2};


%% EXPLORE RAW SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);
plot_NIR(str2num(Xmnir2021.v), Xmnir2021.d) ;
plot_NIR(str2num(Xmnir2022.v), Xmnir2022.d) ;

%% PREPROCESSING - PART I
% -------------------------------------------------------------------------

% Kaiser 2021 -------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xkaiser2021.v)-400));  % 520
[~,i2] = min(abs(str2num(Xkaiser2021.v)-1800)); % 1800
Xkaiser2021 = selectcol(Xkaiser2021,i1:i2);

% Kaiser 2022 -------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xkaiser2022.v)-400));  % 520
[~,i2] = min(abs(str2num(Xkaiser2022.v)-1800)); % 1800
Xkaiser2022 = selectcol(Xkaiser2022,i1:i2);

% MikroNIR ----------------------------------------------------------------
% Keep full region

clear i1 i2

%% EXPLORE NEW SPECTRAL REGIONS
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);
plot_NIR(str2num(Xmnir2021.v), Xmnir2021.d) ;
plot_NIR(str2num(Xmnir2022.v), Xmnir2022.d) ;


%% SPIKE REMOVAL FOR RAMAN SPECTRA 
%  ------------------------------------------------------------------------
% Run to check, but no spikes in these spectra.

th = 20; % threshold for spike detection
[Xkaiser2021, ~] = spikefix_whitaker_multi(Xkaiser2021,2,2,1,th,1); 

th = 20; % threshold for spike detection
[Xkaiser2022, ~] = spikefix_whitaker_multi(Xkaiser2022,2,2,1,th,1); 


%% SMOOTHING RAMAN SPECTRA 
% Savitzky-Golay
% -------------------------------------------------------------------------
wd = 15; % window size
polorder = 2; % polynomial order
derorder =  0; % derivative order 0 = smoothing

% Kaiser 2021 -------------------------------------------------------------
Xkaiser2021 = saisir_derivative(Xkaiser2021,polorder,wd,derorder);
nvar = length(Xkaiser2021.v);
Xkaiser2021 = selectcol(Xkaiser2021,8:(nvar-7));% Remove edge effects (the edge points on each side)

% Kaiser 2022 -------------------------------------------------------------
Xkaiser2022 = saisir_derivative(Xkaiser2022,polorder,wd,derorder);
nvar = length(Xkaiser2022.v);
Xkaiser2022 = selectcol(Xkaiser2022,8:(nvar-7));% Remove edge effects (the edge points on each side)

%% EXPLORE SMOOTHED SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);


%% SEPARATE INSTRUMENTS AND SUBSETS
% -------------------------------------------------------------------------

% Different region versions
X_instrument_subsets = {Xkaiser2021, Xkaiser2022, Xmnir2021, Xmnir2022};
X_subsets = {};
n_subsets = {};

for i = 1:size(X_instrument_subsets,2)

    X = X_instrument_subsets{i};
    % Separate by muscle/part
    [X_sunnyside, ~] = select_from_identifier(X,5,'A'); 
    [X_shadowside, ~] = select_from_identifier(X,5,'B'); 

    % Gather all data subsets ------------------------------------------------- 
    nsets = 3;
    X_subsets(end+1:end+nsets) = {X_sunnyside, X_shadowside, X};
    n_subsets(end+1 : end+nsets , 1) = {'Sunny side'; 'Shadow side'; 'All'};

end



keep X_subsets n_subsets plot_colors References

%% CONNECT REFERENCE DATA WITH RAMAN SPECTRA
% -------------------------------------------------------------------------

ref_idpos = [1:4,9:14]; % string position in reference id tags
x_idpos = [1:4,9:14]; %string position in spectral/data blocks 
reference_names = References.i(:,ref_idpos);
nref = size(References.d,2); % Number of reference analyses
Y_subsets = {};

for i_set = 1:length(X_subsets)  % For each subset, create a reference block with row-to-row correspondence with spectral blocks

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

    % Check if all spectra has a reference value assigned in Y block
    nanflag = any(isnan(Y_subsets{i_set}.d(:))); 
    if nanflag 
        warning(['NAN values were detected in Y_subsets{',num2str(i_set),'}.'])
    end

end



%% BASELINE CORRECTION ONLY (by ALS) FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

X_prep_als_subsets = {};
lambda = 4 ;% Smoothing parameter % 4
p = 0.001;

for i = 1:6 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    X_prep_als_subsets{i} = Subset_prep;
    
    % % Control plot
    % plot_Raman(str2num(Subset_prep.v), Subset_prep.d);

end

%% COMPARE 2021 AND 2022 RAMAN SPECTRA (BEFORE AND AFTER INSTRUMENT SERVICE)
% -------------------------------------------------------------------------
X2021 = X_prep_als_subsets{3};
X2022 = X_prep_als_subsets{6};

ph = plot_Raman(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));

leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)

%% BASELINE CORRECTION (by ALS) and SNV FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

X_prep_alssnv_subsets = {};
lambda = 4; % Smoothing parameter % 4
p = 0.001;

for i = 1:6 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    Subset_prep = saisir_snv(Subset_prep);
    X_prep_alssnv_subsets{i} = Subset_prep;
    
    % % Control plot
    % plot_Raman(str2num(Subset_prep.v), Subset_prep.d);
    % color_by(Y_subsets{i}.d)

end


%% COMPARE 2021 AND 2022 RAMAN SPECTRA (BEFORE AND AFTER INSTRUMENT SERVICE)
% -------------------------------------------------------------------------
X2021 = X_prep_alssnv_subsets{3};
X2022 = X_prep_alssnv_subsets{6};

ph = plot_Raman(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));

% Sugar related peaks:
xline(1458)
xline(1265)
xline(1078)
xline(1124)
xline(918)
xline(629)
xline(493)
xline(834)

leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)



%% SNV FOR NIR SPECTRA
% -------------------------------------------------------------------------

X_prep_snv_subsets = {};

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_snv(Subset); % SNV correction
    X_prep_snv_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    %color_by(Y_subsets{i}.d(:,1));

end


%% OVERALL VALIDATION PERFORMANCES - LOOCV GENERAL MODELS / RAMAN KAISER
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

idpos = [1:3,9:14];
aopt = 20;
target = 1;
[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X_alssnv, Y, target, idpos, aopt,'PlotLimits',[1 20]);

plot_regcoef(X_alssnv.v,squeeze(Regcoeff.d(:,:,4)),0.5, 'findPeaks','no');

%% VALIDATION BY MONTH / ALL COMBINATIONS VALIDATION SCHEME / RAMAN KAISER
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

% Validate from identifier 
aopt = 15;
[FIG, Regcoeff, Performance] = PLSR_allcombos_validation_from_identifier(X_alssnv,Y,1,...
                                9:14,aopt,[1:3,9:14], 'aoptMode','calCV', 'aoptAlg','WM', ...
                                'plotLimits', [1 20] );
lgd = findobj('type','legend');
delete(lgd)

% Plot performance summary
[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 1.7])
set(ah{2}, 'YLim',[-2 2])
set(ah{3}, 'YLim',[0, 1.5])

plot_regcoef(Regcoeff.v, Regcoeff.d, 0.5,'findPeaks','no','plotColors', plot_colors)
xlabel('Raman shift (cm^{-1})', 'FontSize',18)
xline(1458, 'Color', [0.5 0.85 0.5],'Linewidth', 3) % Indicate sugar related peaks
xline(1265, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(1078, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(1124, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(918, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(629, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(493, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(834, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
legend(cellstr(Regcoeff.i), 'Location', 'Eastoutside')


%% OVERALL VALIDATION PERFORMANCES - LOOCV GENERAL MODELS / MICRONIR
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

idpos = [1:3,9:14];
aopt = 20;
target = 1;
[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X_snv, Y, target, idpos, aopt,'PLotLimits',[1 20]);

figs = findobj('type','Figure');
FIG1 = figs(2);
FIG2 = figs(1);

plot_regcoef(X_snv.v,squeeze(Regcoeff.d(:,:,12)),0.5, 'findPeaks','no')

%% VALIDATION BY MONTH / ALL COMBINATIONS VALIDATION SCHEME / MICRONIR
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});


% Validate from identifier 
aopt = 5; % 15
[FIG, Regcoeff, Performance] = PLSR_allcombos_validation_from_identifier(X_snv,Y,1, ...
                                9:14,aopt,[1:3,9:14],'aoptMode','calCV','aoptAlg','WM',...
                                'plotLimits', [1 20]); 
lgd = findobj('type','legend');
delete(lgd)

% Plot performance summary
[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 1.7])
set(ah{2}, 'YLim',[-2 2])
set(ah{3}, 'YLim',[0, 1.5])

plot_regcoef(Regcoeff.v, Regcoeff.d, 0.95,'plotColors',plot_colors,'findpeaks','no')
xline(1200,'LineWidth',30,'Color',[0.5 0.85, 0.5])
xline(1437, 'LineWidth',30,'Color',[0.5 0.85, 0.5])
legend(cellstr(Regcoeff.i),'Location','Eastoutside')
xlabel('Wavelength (nm)','FontSize',18)


%% BOOTSTRAP VALIDATION (AS A FUNCTION OF NUMBER OF SAMPLES) / RAMAN KAISER
% -------------------------------------------------------------------------
% Bootstrap validation from 2021 (cal) to 2022 (val) data, with 
% increasing number of calibration samples. Tha model complexity is chosen 
% in three different ways.

% Merge 2021 and 2022 data blocks
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

% These following snippets are time consuming to run, but could be faster
% by decreasing the variations in number of calibration samples or the
% number of draws.
rng default
nsamples = [3:3:130];
ndraws = 1000;
target = 1; % Target value /col in Y
ncompmax = 5; % 15

[Result_summary_Raman_testrmsepcor, Result_Perdraw_Raman_testrmsepcor] = ...
    PLSR_bootstrap_nsamplesval(X_alssnv,Y,target,ndraws, nsamples, 9:12, ncompmax,'aoptAlgMetric','RMSEP_corr');

[Result_summary_Raman_testrmsep,Result_Perdraw_Raman_testrmsep] = ...
    PLSR_bootstrap_nsamplesval(X_alssnv,Y,target,ndraws, nsamples, 9:12, ncompmax,'aoptAlgMetric','RMSEP'); 

[Result_summary_RamanCV,Result_Perdraw_RamanCV] = ...
    PLSR_bootstrap_nsamplesvalCV(X_alssnv,Y,target,ndraws, nsamples, 9:12,[1:3,9:14],ncompmax);


% Plot mean RMSEPs as a function of number of samples
Result_summary_Raman = Result_summary_RamanCV;


%% BOOTSTRAP VALIDATION (AS A FUNCTION OF NUMBER OF SAMPLES) / MICRONIR
% -------------------------------------------------------------------------
% Bootstrap validation from 2021 (cal) to 2022 (val) data, with 
% increasing number of calibration samples. Tha model complexity is chosen 
% in three different ways.

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12}); 
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

% These following snippets are time consuming to run, but could be faster
% by decreasing the variations in number of calibration samples or the
% number of draws.
rng default
nsamples = [3:3:130];
ndraws = 1000;
target = 1; % Target value /col in Y
ncompmax = 5; % 15

[Result_summary_NIR_testrmsepcor, Result_Perdraw_NIR_testrmsepcor] = ...
    PLSR_bootstrap_nsamplesval(X_snv,Y,target,ndraws, nsamples, 9:12, ncompmax,'aoptAlgMetric', 'RMSEP_corr'); 

[Result_summary_NIR_testrmsep, Result_Perdraw_NIR_testrmsep]  = ... 
    PLSR_bootstrap_nsamplesval(X_snv,Y,target,ndraws, nsamples, 9:12, ncompmax,'aoptAlgMetric', 'RMSEP'); 

[Result_summary_NIRCV, Result_Perdraw_NIR_rmsecv]  = ...
    PLSR_bootstrap_nsamplesvalCV(X_snv,Y,target,ndraws, nsamples, 9:12,[1:3,9:14],ncompmax); 

% Plot mean RMSEPs as a function of number of samples
Result_summary_NIR = Result_summary_NIRCV;


%% COMPARE RAMAN AND NIR BY BOOTSTRAP METHOD, FOR #SAMPLES = 102 
% (Performances stabilized for that large number. and still leaves room
% for variation in drawn samples for the Brix ranges)
% -------------------------------------------------------------------------

row_nsamp102 = 34; % Double check index before running

% Violin plots ------------------------------------------------------------
errortype = 3; % 8,6,4,3 (index in the Results summary struct)
AOPTMAX = input('Enter the maxium component restriction: ');

ViolinBootTestRMSEPcorr = Result_Perdraw_Raman_testrmsepcor.d(34,:,errortype); 
ViolinBootTestRMSEP = Result_Perdraw_Raman_testrmsep.d(34,:,errortype);
ViolinBootRMSECV = Result_Perdraw_RamanCV.d(34,:,errortype);
NIRViolinBootTestRMSEPcorr = Result_Perdraw_NIR_testrmsepcor.d(34,:,errortype);
NIRViolinBootTestRMSEP = Result_Perdraw_NIR_testrmsep.d(34,:,errortype);
NIRViolinBootRMSECV = Result_Perdraw_NIR_rmsecv.d(34,:,errortype);

BootStrapModes = {'Raman_{RMSEPcorr}';'Raman_{RMSEP}';'Raman_{RMSECV}';...
                  'NIR_{RMSEPcorr}';'NIR_{RMSEP}';'NIR_{RMSECV}'};

figure
vs = violinplot([ViolinBootTestRMSEPcorr' ViolinBootTestRMSEP' ViolinBootRMSECV' NIRViolinBootTestRMSEPcorr' NIRViolinBootTestRMSEP' NIRViolinBootRMSECV'],...
                BootStrapModes,'MarkerSize', 1,'MedianMarkerSize', 45, 'QuartileStyle','shadow');
ylabel(Result_Perdraw_NIR_rmsecv.k(errortype,:))

box off 
grid on
set(gcf,'renderer','Painters')




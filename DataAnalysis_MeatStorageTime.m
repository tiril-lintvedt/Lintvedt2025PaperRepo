% =========================================================================
% PAPER // Raman and NIR spectroscopy: A discussion of calibration robustness 
% for food quality measurements through two case studies
%
% Supplemental case: Meat storage time  
% =========================================================================
%
% ------------------------------------------------------------------------
% 
% Toolboxes required: Saisir, GSTools, My_toolbox
%
%
% Author: Tiril Lintvedt
% Nofima, Raw materials and process optimisation
% email address: tiril.lintvedt@nofima.no
% Website: 
% Oct 2023; 
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
              0.2 0.7 0.7    ];

set(0,'defaultfigurecolor',[1 1 1])

addpath(genpath('Data'))
addpath(genpath('Saisir'))
addpath(genpath('GSTools'))
addpath(genpath('My_toolbox'))

%% DATA IMPORT
% -------------------------------------------------------------------------


% MM 
spec = GSImportspec('Data\Meat\20240305 Meat oxygenation supplementary',0);
[X, Xmeta] = extractSPC(spec); % Spectra and Meta data



%% PREPROCESSING - PART I
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(X.v)-515));  % 520
[~,i2] = min(abs(str2num(X.v)-1870)); % 1800 % 1770
X = selectcol(X,i1:i2);               


% Remove non-meat spectra
X_cyclo = select_from_identifier(X,1,'Cyclo');
X = delete_from_identifier(X,1,'Cyclo');




%% EXPLORE NEW SPECTRAL REGIONS
% -------------------------------------------------------------------------

plot_Raman(str2num(X.v), X.d);


%% SPIKE REMOVAL FOR RAMAN SPECTRA 
%  ------------------------------------------------------------------------

% Remove damaged pixels/spikes MM
dpix =[660:664, 950:954, 1140:1144, 1584:1588, 1768:1800, 1980:1983, 1996:1999]; 
[~,di] = min(abs(str2num(X.v)-dpix)); % 1800
X = deletecol(X,di);
th = 6.5; % threshold for spike detection
[X, ~] = spikefix_whitaker_multi(X,2,1,th,1); % double check this later



%% EXPLORE SPIKE FREE SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(X.v), X.d);


%% SMOOTHING RAMAN SPECTRA
% Savitzky-Golay
% -------------------------------------------------------------------------
wd = 9; % window size
polorder = 2; % polynomial order
derorder =  0; % derivative order 0 = smoothing

X = saisir_derivative(X,polorder,wd,derorder); 
nvar = length(X.v);
X = selectcol(X,5:(nvar-4)); % Remove edge effects (the edge points on each side)

%% EXPLORE SMOOTHED SPECTRA 
% -------------------------------------------------------------------------

hour = str2num(char(strrep(cellstr(X.i(:,18:20)),'e','')));
plot_Raman(str2num(X.v), X.d);
ch = color_by(hour);
title(ch,'Hours')
grid on


%  Look only at day 1
minutes = 60.*hour + str2num(char(cellstr(X.i(:,22:23)))); 
Xday1  = select_from_identifier(X,19,'0');
hour = str2num(char(strrep(cellstr(Xday1.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(Xday1.i(:,22:23)))); 
plot_Raman(str2num(Xday1.v), Xday1.d);
ch = color_by(minutes);
title(ch,'Minutes')
grid on


%% SEPARATE INSTRUMENTS AND SUBSETS
% -------------------------------------------------------------------------

% Different region versions
X_instrument_subsets = {X};
X_subsets = {};
n_subsets = {};

for i = 1:size(X_instrument_subsets,2)
    
    X = X_instrument_subsets{i};

    % Separate by muscle/part
    [X_Bankekjott, ~] = select_from_identifier(X,1,'B'); 



    % Gather all data subsets ------------------------------------------------- 
    nsets = 1;
    X_subsets(end+1:end+nsets) = {X_Bankekjott};
    n_subsets(end+1 : end+nsets , 1) = {'Bankekjøtt (intakt muskel)'};

end



%% BASELINE CORRECTION ONLY (by ALS) FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

X_prep_als_subsets = {};
Fbaseline_subsets = {};
lambda = 3.5; % Smoothing parameter % 4
p = 0.001;

for i = 1 
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    X_prep_als_subsets{i} = Subset_prep;

    Baselines.i = Subset.i;
    Baselines.v = Subset.v;
    Baselines.d = baseline;
    Fbaseline_subsets{i} = Baselines;
    
    % % Control plot
    % plot_Raman(str2num(Subset_prep.v), Subset_prep.d);

end



%% EXPLORE BASELINE CORRECTED SPECTRA OVER TIME
% -------------------------------------------------------------------------

% Color by fluorescence levels, all storage days
isubset = 1;
plot_Raman(str2num(X_prep_als_subsets{isubset}.v), X_prep_als_subsets{isubset}.d); 
ch = color_by(mean(Fbaseline_subsets{isubset}.d,2));
title(ch,{'Fluorescence';'ind'})
grid on

% Color by fluorescence levels, storage day 1
Xday1  = select_from_identifier(X_prep_als_subsets{isubset},19,'0');
Fday1 = select_from_identifier(Fbaseline_subsets{isubset},19,'0');
plot_Raman(str2num(Xday1.v), Xday1.d) 
ch = color_by(mean(Fday1.d,2));
title(ch,{'Fluorescence ind.'})
grid on


% Color by time
hour = str2num(char(strrep(cellstr(X.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(X.i(:,22:23)))); 

% All days
isubset = 1;
plot_Raman(str2num(X_prep_als_subsets{isubset}.v), X_prep_als_subsets{isubset}.d) ;
color_by(minutes)


%% BASELINE CORRECTION ONLY (MODPOLY) FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

X_prep_modpol_subsets = {};
Fbaseline_modpol_subsets = {};

for i = 1 
    Subset = X_subsets{i};
    [Subset_prep, Baselines] = saisir_modpoly(Subset, 'order', 10); % correct baseline
    X_prep_modpol_subsets{i} = Subset_prep;
    Fbaseline_modpol_subsets{i} = Baselines;
    
    % % Control plot
    % plot_Raman(str2num(Subset_prep.v), Subset_prep.d);
end



%% EXPLORE BASELINE CORRECTED SPECTRA 
% -------------------------------------------------------------------------


% Color by fluorescence levels, all storage days
isubset = 1;
plot_Raman(str2num(X_prep_modpol_subsets{isubset}.v), X_prep_modpol_subsets{isubset}.d) ;
color_by(mean(Fbaseline_modpol_subsets{isubset}.d,2));

% Color by fluorescence levels, storage day 1
Xday1  = select_from_identifier(X_prep_modpol_subsets{isubset},19,'0');
Fday1 = select_from_identifier(Fbaseline_modpol_subsets{isubset},19,'0');
plot_Raman(str2num(Xday1.v), Xday1.d); 
ch = color_by(mean(Fday1.d,2));
title(ch,{'Fluorescence ind.'})
grid on


% Color by time
hour = str2num(char(strrep(cellstr(X.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(X.i(:,22:23)))); 

% All days
isubset = 1;
plot_Raman(str2num(X_prep_modpol_subsets{isubset}.v), X_prep_modpol_subsets{isubset}.d);
ch = color_by(minutes);
title(ch,{'Minutes'})
grid on


%% FLUORESCENCE LEVEL AS A FUNCTION OF TIME 
% -------------------------------------------------------------------------
F_baseline = Fbaseline_subsets{1};
F_baseline = select_from_identifier(F_baseline,19,'0'); % Day1 only

hour = str2num(char(strrep(cellstr(F_baseline.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(F_baseline.i(:,22:23)))); % all days
unique_id = char(unique(cellstr(F_baseline.i(:,12:13)))); % Identify all unique id tags in position colour_pos

phandles = {};
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
hold on
for  i = 1:length(unique_id)
    [Fbslsub, ind] = select_from_identifier(F_baseline,12, unique_id(i,:));
    Fbslavgsub = saisir_mean(Fbslsub);
    ph = plot(minutes(ind),Fbslavgsub.d,'o','Color', plot_colors(i,:));
    phandles{i} = ph;
end

set(gcf,'Color',[1 1 1])
box off
leg = legend([phandles{1},phandles{2}, phandles{3}], unique_id);
title(leg, 'Sample', 'Fontsize',14)
xlim([-10 600])
ylabel('Fluorescence indicator','FontSize',18)
xlabel('Minutes','FontSize',18)
grid on

%% FLUORESCENCE LEVEL AS A FUNCTION OF TIME 
% -------------------------------------------------------------------------
F_baseline = Fbaseline_subsets{1};
%F_baseline = select_from_identifier(F_baseline,19,'0'); % DAY 1-8

hour = str2num(char(strrep(cellstr(F_baseline.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(F_baseline.i(:,22:23)))); % all days
unique_id = char(unique(cellstr(F_baseline.i(:,12:13)))); % Identify all unique id tags in position colour_pos

phandles = {};
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/2 scrsz(4)/2.5])  
hold on
for  i = 1:length(unique_id)
    [Fbslsub, ind] = select_from_identifier(F_baseline,12, unique_id(i,:));
    Fbslavgsub = saisir_mean(Fbslsub);
    ph = plot(hour(ind),Fbslavgsub.d,'o','Color', plot_colors(i,:));
    phandles{i} = ph;
end

set(gcf,'Color',[1 1 1])
box off
leg = legend([phandles{1},phandles{2}, phandles{3}], unique_id);
title(leg, 'Sample', 'Fontsize',14)
xlim([-3 200])
ylabel('Fluorescence indicator','FontSize',18)
xlabel('Hours','FontSize',18)
grid on



%% FAT PREDICTION FOR DIFFERENT TIME GROUPS
% -------------------------------------------------------------------------

% Load dummy fat model from Meat case study
load('Data\TotFatModelHomMM.mat')

B = TotFatModel.B;
B = B(2:end); % Random, but need correspondance with reg coef - only initial check.
b0 = TotFatModel.b0;

% day 1, only:
X  = select_from_identifier(X_prep_als_subsets{1},[19],'0');
F = select_from_identifier(Fbaseline_subsets{1},19,'0');
F.d = mean(F.d,2);

% Sample 1 only
samplid = '1';
X  = select_from_identifier(X,13,samplid);
F = select_from_identifier(F,13,samplid);

% Predictions
Ypred.d = b0 + X.d*B;
Ypred.i = X.i;

% Color by time
hour = str2num(char(strrep(cellstr(X.i(:,18:20)),'e','')));
minutes = 60.*hour + str2num(char(cellstr(X.i(:,22:23)))); 

figure; hold on;
for i = 1:length(Ypred.d)
    t = minutes(i);
    plot(t,Ypred.d(i),'o')
end

xlabel('Minutes','FontSize',18)
ylabel('Predicted fat %','FontSize',18)
cb = color_by(F.d);
title(cb, 'Fluorescence level')


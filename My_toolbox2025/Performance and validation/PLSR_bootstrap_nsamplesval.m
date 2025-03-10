function [Result_summary,MetricsPerDraw] = PLSR_bootstrap_nsamplesval(X, Y, target, ndraws, nsamples, idpos, ncomp, varargin)
% ------- Bootstrap validation as a function of number of samples ------s---
% A bootstrap validation version for PLSR validation, with separate cal and
% test sets, where model complexity choice is based on performance in the
% test set. The validation is done for different variations in number of
% calibration samples, ensuring a balanced design in each sample selection
% (based on three groupings of the target value).
% -------------------------------------------------------------------------
%
%   INPUT:
%               X         -  Saisir data structure
%               Y         -  Saisir data structure of reference values 
%               target    -  Column number in block Y to establish models 
%                            for (only one!)
%               ndraws    - Number of random draws for picking calibration
%                           samples
%               nsamples  - Vecor, list of number of samples to look at. The 
%                           number in the end might not be exact. Depends on
%                           divisibility with 3 (Low, Med, High vals).
%               idpos     - Id position (indices) in strings used to define
%                           the main calibration and test set.
%               ncomp      - Predefined number of PLS components to use
%   
%   OUTPUT:
%               Result_summary  -   Summary of all important results 
%               MetricsPerDraw  -  All performance metrics across the
%                                  variations in number of sample included 
%                                  and corresponding draw numbers 
%
% -------------------------------------------------------------------------
% EXAMPLE CALL
% -------------------------------------------------------------------------
% Result_summary = PLSR_montecarlo_nsamplesval(X_instrument_subsets{1}, Y{1},...
%                                               1, 10, [10 20 30 40 50], 13:14, 5)
% * ncomp is an optional parameter, but must be present if 'aoptAlg' is set
%   to 'STATIC' 
% * aoptAlg : 'STATIC' / 'WM' (default) / 'ALT2'
%
% -------------------------------------------------------------------------
% REFERENCES
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SETTINGS 
% -------------------------------------------------------------------------

% Main plot colors
plot_colors = [0 0.5 0.9   ; 
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

plot_markers  = {'^','o','square','diamond','*','x','hexagram','.',...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.'};

set(0,'defaultfigurecolor',[1 1 1])

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultNcomp = 10;
defaultAoptAlg = 'WM'; % // 'ALT2', 'STATIC'
defaultMetric = 'RMSEP'; 
%defaultPlotLimit = [];
defaultSets = [];

p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   validNumLen2= @(x) isnumeric(x) && (length(x)==2) ;
   expectedAlgs = {'STATIC', 'WM', 'ALT2'};
   expectedAlgMetrics = {'RMSEP','RMSEP_corr'};
   
   addRequired(p,'X'); % positional arg
   addRequired(p,'Y'); % positional arg
   addRequired(p,'target',validScalarPosNum); % positional arg
   addRequired(p,'ndraws'); % positional arg
   addRequired(p,'nsamples'); % positional arg
   addRequired(p,'idpos'); % positional arg
   addOptional(p,'ncomp',defaultNcomp,validScalarPosNum);
   %addParameter(p,'plotLimits',defaultPlotLimit,validNumLen2);
   addParameter(p,'aoptAlg',defaultAoptAlg,  @(x) any(validatestring(x,expectedAlgs))); % Name Value pair
   addParameter(p,'aoptAlgMetric',defaultMetric,  @(x) any(validatestring(x,expectedAlgMetrics))); % Name Value pair
   addParameter(p,'Cal_Test',defaultSets, @(x) ischar(x)); % Name Value pair % Id position names to control whats cal and test set respectively


   parse(p,X, Y, target, ndraws, nsamples, idpos, ncomp, varargin{:});
   
 aoptmetric = p.Results.aoptAlgMetric;

 if ~isempty(p.Results.Cal_Test)
     valsets = split(cellstr(p.Results.Cal_Test),'_');
     calset = char(valsets{1});
     testset = char(valsets{2});

 else
    unique_id = char(unique(cellstr(X.i(:,idpos)))); % Unique id
    if size(unique_id,1) ~= 2; error('The number of unique groups must be 2 to define cal and test set.'); end
    calset = unique_id(1,:);
    testset = unique_id(2,:);
 end

% -------------------------------------------------------------------------
% DEFINE CAL. AND VAL. DATA SETS ------------------------------------------
% -------------------------------------------------------------------------

% Divide based on id tag position
Xcal = select_from_identifier(X,idpos,calset);
Ycal = select_from_identifier(Y,idpos,calset);
Xtest = select_from_identifier(X,idpos,testset);
Ytest = select_from_identifier(Y,idpos,testset);

% Define high/medium/low value groups based on Reference values, for the
% balanced design for the calibration sample draw
[Xcal_segments, Ycal_segments, nsampseg] = group_from_yvalue(Xcal,Ycal,1,3);

% check that num samples in each segment is larger than max number of
% samples to draw.  
if min(nsampseg) < max(nsamples./3)
    warning('Number of samples in one or more of the value ranges are too small. ')
    disp(['Num samples in value groups: ', num2str(nsampseg)])
    disp(['Num samples to draw in each value group (for increasnig sample number): ', num2str(ceil(nsamples./3))])
end

% -------------------------------------------------------------------------
% BALANCED RANDOM DRAW, MODELLING AND VALIDATION 
% -------------------------------------------------------------------------

nsampmetrics_descrStats.d = zeros(7,8,length(nsamples));
nsampmetrics_descrStats.i = ['Mean           ';'Median         '; 'Max            ';...
                             'Min            ';'Std            ';'5th percentile ';...
                             '95th percentile']; % descriptive statistics
nsampmetrics_descrStats.v = ['RMSEP        ';'R2P          '; 'BIAS         ';
                             'SLOPE        ';'INTERCEPT    '; 'LV           ';
                             'R2P CORR     ';'RMSEP CORR   '];
nsampmetrics_descrStats.z = [repelem('Num samples: ',length(nsamples),1), num2str(nsamples')];

ci = 0;
for i = 1:length(nsamples)
    ci = ci +1;
    nsamp = nsamples(i);

    ndrawsmetrics.d = zeros(ndraws,8);
    ndrawsmetrics.i = [repelem('draw ',ndraws,1), num2str((1:ndraws)')];
    ndrawsmetrics.v = ['RMSEP        ';'R2P          '; 'BIAS         ';
                       'SLOPE        ';'INTERCEPT    '; 'LV           ';
                       'R2P CORR     ';'RMSEP CORR   '];
    cj = 0;
    for j = 1:ndraws        
        cj = cj+1;

        % Sample selection ------------------------------------------------

        % Random sample selection from each ref. value range group
        Xcal_rand_subsets = cell(size(Xcal_segments));
        Ycal_rand_subsets = cell(size(Ycal_segments));

        for k = 1:length(Xcal_segments)
            nsubsamp = size(Xcal_segments{k}.d,1);
            randind = randi(nsubsamp,ceil(nsamp/3),1); % generate random sample indices within the right index range 
            Xcal_rand_subsets{k} = selectrow(Xcal_segments{k},randind);
            Ycal_rand_subsets{k} = selectrow(Ycal_segments{k},randind);
        end

        % Merge the three sub value ranges
        Xcal_randsel = saisir_rowmerge(Xcal_rand_subsets{:});
        Ycal_randsel= saisir_rowmerge(Ycal_rand_subsets{:});
        
        % Calibration -----------------------------------------------------
       
        % Make PLSR calibration model
        [b0,B,T,W,P,q] = pls_nipals(Xcal_randsel.d,Ycal_randsel.d(:,target),ncomp); 
        
        % Validation ------------------------------------------------------

        % Predict test set samples
        YpredTest = b0 + Xtest.d*B;

        % Calculate RMSE for all number of comps, incl zero comp
        RMSEP = cal_rmse(Ytest.d(:,target), YpredTest);                 
        ypred_zerocomp = mean(Ycal_randsel.d(:,target)); 
        rmse_zerocomp = cal_rmse(repelem(ypred_zerocomp,length(Ytest.d(:,target)))', Ytest.d(:,target));
        R2P =  cal_r2(Ytest.d(:,target), YpredTest);     
        BIAS = cal_bias(Ytest.d(:,target), YpredTest);
        slope_intercept = cal_slope(Ytest.d(:,target), YpredTest) ;

        R2P_corr =  R2_corr(Ytest.d(:,target), YpredTest);
        RMSEP_corr = rmse_corr(Ytest.d(:,target), YpredTest);

        % Find optimal number of components from test set by given strategy
        if  strcmp(p.Results.aoptAlg, 'WM')
            aopt = find_aopt(eval(aoptmetric)); % Based RMSEP/RMSEP_corr, Westad&Martens 
        elseif strcmp(p.Results.aoptAlg, 'ALT2')
            aopt = find_aopt([rmse_zerocomp RMSEP], 0.02, 2);
            aopt = p.Results.ncomp;
        end
        
        % Store metrics ---------------------------------------------------

        % performance metrics for optimal number of components 
        ndrawsmetrics.d(j,:) = [RMSEP(:,aopt), R2P(:,aopt), BIAS(:,aopt), ...
                                slope_intercept(1,aopt), slope_intercept(2,aopt), ...
                                aopt, R2P_corr(:,aopt), RMSEP_corr(:,aopt)];


        % Store metrics for each draw (i:ndraws)
        MetricsPerDraw.d(ci,cj,:) = [RMSEP(:,aopt), R2P(:,aopt), BIAS(:,aopt), ...
                                    slope_intercept(1,aopt), slope_intercept(2,aopt), ...
                                    aopt, R2P_corr(:,aopt), RMSEP_corr(:,aopt)];

        % -----------------------------------------------------------------

    end

    MetricsPerDraw.i = [repelem('Number of samples :', length(nsamples'),1) num2str(nsamples')];
    MetricsPerDraw.j = [repelem('Draw number : ',length(1:ndraws),1) num2str((1:ndraws)')];
    MetricsPerDraw.k = ['RMSEP(AOPT)        ';'R2P(AOPT)          '; 'BIAS(AOPT)         ';
                       'SLOPE(AOPT)        ';'INTERCEPT(AOPT)    '; 'LV(AOPT)           ';
                       'R2P CORR(AOPT)     ';'RMSEP CORR(AOPT)   '];    

    % RESULTS DISTRIBUTION ------------------------------------------------
    % ---------------------------------------------------------------------
    
    % Descriptive statistics for the distribution of performance across ndraws
    mean_metrics = mean(ndrawsmetrics.d,1);
    median_metrics = median(ndrawsmetrics.d,1);
    max_metrics = max(ndrawsmetrics.d,[],1);
    min_metrics = min(ndrawsmetrics.d,[],1);
    std_metrics = std(ndrawsmetrics.d,1);
    prct5th_metrics = prctile(ndrawsmetrics.d,5);
    prct95th_metrics = prctile(ndrawsmetrics.d,95);
    
    % Merge in one matrix
    DescrStats_i = [mean_metrics;median_metrics;max_metrics;min_metrics;...
                    std_metrics;prct5th_metrics;prct95th_metrics];

    nsampmetrics_descrStats.d(:,:,i) = DescrStats_i;

end


 Result_summary.nsamples = nsamples;
 Result_summary.ndraws = ndraws;
 Result_summary.nsamplesPerValueGroup = nsampseg; 
 Result_summary.nsamples2drawPerValueGroup = ceil(nsamples./3); 
 Result_summary.DescrStatsMetrics = nsampmetrics_descrStats;
 

% PLOTS -------------------------------------------------------------------
% -------------------------------------------------------------------------

% Plot mean RMSEPs as a function of number of samples
figure
hold on
plot(nsamples, squeeze(Result_summary.DescrStatsMetrics.d(2,1,:)),'-o')
plot(nsamples, squeeze(Result_summary.DescrStatsMetrics.d(6,1,:)),'--')
plot(nsamples, squeeze(Result_summary.DescrStatsMetrics.d(7,1,:)),'--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median RMSEP', 'FontSize',18)
box off

end
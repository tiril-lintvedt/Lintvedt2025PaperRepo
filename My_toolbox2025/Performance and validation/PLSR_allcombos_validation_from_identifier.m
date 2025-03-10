function [FIG, Regcoeff, Performance, Residuals] = PLSR_allcombos_validation_from_identifier(X, Y, target, idpos, ncomp,reps_idpos, varargin)
% ---------------------- All-combinations-validation ----------------------
% A version of PLSR validation where all unique groups defined by the saisir 
% id string are used as a calibration set and test set, testing all possible 
% combinations. This is different from normal crossvalidation because each 
% sample is predicted several times. Finds the optimal number of components 
% from optimal test set performance. 
% 
% NB! Works only for one target reference value at a time.
% -------------------------------------------------------------------------
%
%   INPUT:
%               X           -   Saisir data structure
%               Y           -   Saisir data structure of reference values 
%               target      -   Column number in block Y to establish models 
%                               for (only one!)
%               idpos       -   Id position (indices) in strings used to 
%                               define validation segments     
%               reps_idpos  -   Id position (indices) in strings used to 
%                               define what measurements belong to same sample.
%               aopt        -   Predefined number of PLS components to use
% 
%   OUTPUT:
%               FIG         - Figure handle for predicted vs target plot
%               Regcoeff    - Regression coefficients per segment
%               Performance - Struct, summarry of performance metrics
%               Residuals   - Prediction residuals
% -------------------------------------------------------------------------
% EXAMPLE CALL
% -------------------------------------------------------------------------
%
% [FIG, Regcoeff, Performance] = PLSR_validate_from_identifier2(X, Y, target, ...
%                                                               idpos, ncomp,'aoptAlg', 'WM')
%  
% * ncomp is an optional parameter, but must be present if 'aoptAlg' is set
%   to 'STATIC' 
% * aoptAlg : 'STATIC' / 'WM' (default) / 'ALT2'
% * reps_idpos is an optional parameter but must be present if 'aoptMode'
%   is set to 'calCV'
% -------------------------------------------------------------------------
% SETTINGS 
% -------------------------------------------------------------------------

% Main plot colors
plot_colors =[0.1 0.1 0.5 ; 
              0.5 0.5 0.5 ;
              0.9 0.5 0   ; 
              0.2 0.7 0.7 ;
              0.6 0.2 0.9 ;
              0.3 0.8 0.3 ;
              0.9 0.3 0.6 ; 
              0   0   0   ;
              0 0.5 0.9   ; 
              0.2 0.5 0.2 ; 
              0.7 0.7 0.2 ; 
              0.9 0.3 0.1 ;
];

% plot_colors = ones(12,3)-linspecer(12,'qualitative') + 0.3*ones(12,3);

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
defaultPlotLimit = [];
defaultAoptMode = 'test';
defaultRepIdPos = [];

p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   validNumLen2= @(x) isnumeric(x) && (length(x)==2) ;
   expectedModes = {'calCV','test'}; 
   expectedAlgs = {'STATIC', 'WM', 'ALT2'};
   expectedAlgMetrics = {'RMSEP','RMSEP_corr'};
   
   addRequired(p,'X'); % positional arg
   addRequired(p,'Y'); % positional arg
   addRequired(p,'target',validScalarPosNum); % positional arg
   addRequired(p,'idpos'); % positional arg
   addOptional(p,'ncomp',defaultNcomp,validScalarPosNum);
   addOptional(p,'reps_idpos',defaultRepIdPos);
   addParameter(p,'plotLimits',defaultPlotLimit,validNumLen2);
   addParameter(p,'aoptAlg',defaultAoptAlg,  @(x) any(validatestring(x,expectedAlgs))); % Name Value pair
   addParameter(p,'aoptMode',defaultAoptMode,  @(x) any(validatestring(x,expectedModes)));
   addParameter(p,'aoptAlgMetric',defaultMetric,  @(x) any(validatestring(x,expectedAlgMetrics))); % Name Value pair
   
   parse(p,X, Y, target, idpos, ncomp, reps_idpos, varargin{:});
   
 aopt = p.Results.ncomp; % Need to add "if exist"
 aoptmetric = p.Results.aoptAlgMetric;
 aoptalg =  p.Results.aoptAlg;
 aoptmode = p.Results.aoptMode;
 plotlim = p.Results.plotLimits;

% -------------------------------------------------------------------------

% Model frames
[m,n] = size(X.d);
mc = aopt ; % NB! Hve you coded parsing for aopt input?a

% Storage for predictions
Ypred.d = zeros([m length(target)]);
Ypred.i = Y.i;

% Determine number of CV segments from identifier and initialize 
unique_id = char(unique(cellstr(X.i(:,idpos)))); % Unique
cv_seg = group_from_identifier(X, idpos); % Assign each spectrum to a CV segment by unique id tag position
unique_group_num = unique(cv_seg); % Unique cv segment numbers
ncvseg = length(unique_group_num); 
Ypred.v = cell(m, ncvseg); % Might be a better way of storing this in the future

Regcoeff.d = zeros(ncvseg,n); 
Regcoeff.v = X.v;
Regcoeff.i = '                          ';

%Residuals.d = zeros(m, ncvseg); 
Residuals.v = '                          '; 

% Initialize predicted versus target figure
scrsz = get(0,'ScreenSize');
FIG = figure('Position',[50 50 scrsz(3)/2.5 scrsz(4)/1.9]);
hold on
phandles = {};
ci = 1; % plotnr counter
group_names = {};

% Cross validation based on id tags
for k = unique_group_num'

    % Sample indices going into calibration set
    indin = (k == cv_seg); % For contigous blocks: setdiff(1:m,[k, k+1]);

    % PLS calibration
    [b0,B,T,W,P,q] = pls_nipals(X.d(indin, :),Y.d(indin,target),mc); 
    
    if strcmp(aoptmode,'calCV')
        % Find aopt
        % Cross validation on calibration set, choose number of comp.
        Xin = selectrow(X,indin);
        Yin = selectrow(Y,indin);
        [~, ~, PerformanceMetrics] = PLSR_CV_from_identifier(Xin, Yin, target,reps_idpos, mc,'showPlot','no','aoptAlgMetric','RMSECV', 'aoptAlg',p.Results.aoptAlg);
        aopt = PerformanceMetrics.lv;
    end

    % Predict cal set from cal model (for all number of comps) and calculate 
    YpredCal.d(indin,:) = b0 + X.d(indin,:)*B;
    YpredCal.i(indin,:) = X.i(indin,:);
    YpredCal.v = ['Number of components (1,2,3, ...)'];
    YRefCal = Y.d(indin,target); % Reference values for calibration set
    
    % Performance metrics cal set 
    RMSEC = cal_rmse(YRefCal, YpredCal.d(indin,:));  
    R2C = cal_r2(YRefCal, YpredCal.d(indin,:));

    % PLS prediction 
    % Use samples from one category to build model, and apply on the other 
    % categories, defined by identifier string positions.
    
    % Keep control of what group was predicted from what other group.
    segsout = find(unique_group_num ~=k); % Determine what are the other group numbers
   
    for j = segsout'
        indout = (j == cv_seg);

        % Prediction for all number of PLS comps in model   
        YpredTest.d(indout,:) = b0 + X.d(indout,:)*B;
        YpredTest.i(indout,:) = X.i(indout,:);
        YpredTest.v = 'Number of components (1,2,3, ...)';
        YRefTest = Y.d(indout,target); % Reference values for calibration
    
        % Calculate RMSE for all number of comps, incl zero comp
        RMSEP = cal_rmse(YRefTest, YpredTest.d(indout,:));          
        RMSEP_corr = rmse_corr(YRefTest, YpredTest.d(indout,:));
        ypred_zerocomp = mean(YRefCal); % Previously this has been done for average of test set
        rmse_zerocomp = cal_rmse(repelem(ypred_zerocomp,length(YRefTest))', YRefTest);
        
        % Find optimal number of components from given strategy
        if strcmp(aoptmode,'test')
            if  strcmp(p.Results.aoptAlg, 'WM')
                aopt = find_aopt(eval(aoptmetric)); % Based RMSEP/RMSEP_corr, Westad&Martens         
            elseif strcmp(p.Results.aoptAlg, 'ALT2')
                aopt = find_aopt([rmse_zerocomp RMSEP], 0.02, 2); % Optimal number of components, Unscr                     
            elseif strcmp(p.Results.aoptAlg, 'STATIC')
                aopt = p.Results.ncomp;
            end
        end
             

        % Optimal model ---------------------------------------------------
        b0_aopt = b0(aopt);
        b_aopt = B(:,aopt);
            
        % Store final regression coefficients (Now only assumes 1 per idtype)        
        descr = [char(unique(cellstr(X.i(indin,idpos)))) '->' char(unique(cellstr(X.i(indout,idpos)))) ' ( LV = ' num2str(aopt) ' )']; % Description of segments used
        Regcoeff.d(ci,:) = b_aopt;
        Regcoeff.i(ci,1:length(descr)) = descr;  %[char(unique(cellstr(X.i(indin,idpos))))];

        % Store predictions and metrics for the OPTIMAL comp number
        Ypred.d(indout,k) = b0_aopt + X.d(indout,:)*b_aopt;
        Ypred.v(indout,k) = cellstr(repelem(descr,sum(indout),1));
               

        % Predicted versus target plot for the OPTIMAL comp number
        p_i = plot(Y.d(indout,target), Ypred.d(indout,k),'.', 'color',plot_colors(ci,:),'Marker',plot_markers{ci},'MarkerSize',5);
        phandles{ci} = p_i;
        plot(-15:50, -15:50, 'color', 'r')
      

        % Calculate performance metrics for OPTIMAL comp number 
        y_test = Y.d(indout,target); % Reference values for current testset
        y_pred = Ypred.d(indout,k); % Predicted values for current test set wit opt. num. comps.
        
        Residuals.d(indout,ci) = y_test - y_pred;
        Residuals.v(ci,1:length(descr)) = char(descr);
        Residuals.i(indout,:) = X.i(indout,:);

        % Residual variance, using zero components 
        ypred_zerocomp = mean(y_test);
        rmse_zerocomp = cal_rmse(repelem(ypred_zerocomp,length(y_test))', y_test);
        
        % Performance metrics 
        RMSEP = cal_rmse(y_test, y_pred); 
        R2P = cal_r2(y_test, y_pred);
        BIAS = cal_bias(y_test, y_pred);
        slope_intercept = cal_slope(y_test, y_pred) ;
        LV = aopt;

        % Performance after slope and bias correction
        RES_corr = residuals_corr(y_test, y_pred);
        RMSEP_corr = rmse_corr(y_test, y_pred);
        R2P_corr = R2_corr(y_test, y_pred);

        % Residuals when bias is corrected
        Residuals.dcor(indout,ci) = RES_corr;

        % Store all final metrics -----------------------------------------
        % INFO: Metrics are store row wise according to nested for-loop 
        Performance.rmsep(k,j) = RMSEP; 
        Performance.rmsepcorr(k,j) = RMSEP_corr; 
        Performance.r2p(k,j) = R2P;
        Performance.r2pcorr(k,j) = R2P_corr; 
        Performance.bias(k,j) = BIAS;
        Performance.slope(k,j) = slope_intercept(1); 
        Performance.lv(k,j) = LV;
        Performance.i(k,:) = [char(unique(cellstr(X.i(indin,idpos)))) '(Calibration)' '->' ];
        Performance.v(j,:) = ['->' char(unique(cellstr(X.i(indout,idpos)))) '(Validation)'];
        

        group_names{ci} = [descr(1:end)];   
        ci = ci+1; 
    end
end

xlim(plotlim) 
ylim(plotlim) 
xlabel('Target','FontSize',18)
ylabel('Predicted','FontSize',18)
leg = legend([phandles{:}],group_names{:},'FontSize',9); % , 'Location', 'Southeast');
title(leg, 'Calibration->Test','FontSize',11)



end
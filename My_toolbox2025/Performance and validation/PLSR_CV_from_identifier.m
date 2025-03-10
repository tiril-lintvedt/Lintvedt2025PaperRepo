function [Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X, Y, target, idpos, ncomp, varargin)
% ----------------- Leave-samples-out  cross validation -------------------
%
% CV with segments defined by samples for identical id strings in given 
% position idpos.
% -------------------------------------------------------------------------
%
%   INPUT:
%               X       -   Saisir data structure, spectra
%               Y       -   Saisir data structure, reference values 
%               target  -   Column number in matrix Y to establish models for (only one)
%               idpos   -   Position (indices) in strings used to define validation segments       
%               aopt    -   Predefined number of PLS components to use
% 
%   OUTPUT:
%               Predictions -   Predictions
%               Regcoeff    -   saisir, Jackknifed regression coefficients for each segment in the CV
%               Performance -   saisir, summary of all performance metric
%
% -------------------------------------------------------------------------
% EXAMPLE CALL
% -------------------------------------------------------------------------
%
% [FIG, Regcoeff, Performance] = PLSR_validate_from_identifier2(X, Y, target, ...
%                                                               idpos, ncomp,'aoptAlg', 'UNSC')
%  
% * ncomp is an optional parameter, but must be present if 'aoptAlg' is set
%   to 'STATIC' 
% * aoptAlg : 'STATIC' / 'WM' (Westad&Martens approach, default) / 'ALT2'
% -------------------------------------------------------------------------
% SETTINGS 
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

plot_markers  = {'^','o','square','diamond','*','x','hexagram','.',...
                 '^','o','square','diamond','*','x','hexagram','.', ...
                 '^','o','square','diamond','*','x','hexagram','.'};

set(0,'defaultfigurecolor',[1 1 1])

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

defaultNcomp = 10;
defaultAoptAlg = 'WM'; % // 'ALT2', 'STATIC'
defaultMetric = 'RMSECV'; 
defaultShowPlot = 'yes';
defaultPlotLimit = [];


p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   validNumLen2= @(x) isnumeric(x) && (length(x)==2) ;
   expectedAlgs = {'STATIC', 'WM', 'ALT2'};
   expectedAlgMetrics = {'RMSECV','RMSECV_corr'};
   expectedShowPlot = {'yes','no'};
   
   addRequired(p,'X'); % positional arg
   addRequired(p,'Y'); % positional arg
   addRequired(p,'target',validScalarPosNum); % positional arg
   addRequired(p,'idpos'); % positional arg
   addOptional(p,'ncomp',defaultNcomp,validScalarPosNum);
   addParameter(p,'showPlot',defaultShowPlot,@(x) any(validatestring(x,expectedShowPlot)));
   addParameter(p,'plotLimits',defaultPlotLimit,validNumLen2);
   addParameter(p,'aoptAlg',defaultAoptAlg,  @(x) any(validatestring(x,expectedAlgs))); % Name Value pair
   addParameter(p,'aoptAlgMetric',defaultMetric,  @(x) any(validatestring(x,expectedAlgMetrics))); % Name Value pair
   
   parse(p,X, Y, target, idpos, ncomp, varargin{:});
   
 aopt = p.Results.ncomp; 
 aoptmetric = p.Results.aoptAlgMetric;
 plotlim = p.Results.plotLimits;
 showPlot = p.Results.showPlot;

% -------------------------------------------------------------------------

% Model frames
[m,n] = size(X.d);
ncomp = min((m-2),aopt) ; 

% Storage for test set predictions (for each number of model comps included)
YpredTest.d = zeros([m ncomp]);
YpredTest.i = Y.i;
YpredTest.v = num2str([1:ncomp]');

% Determine number of CV segments from identifier and initialize 
unique_id = char(unique(cellstr(X.i(:,idpos)))); % Unique
cv_seg = group_from_identifier(X, idpos); % Assign each spectrum to a CV segment by unique id tag position
unique_group_num = unique(cv_seg); % Unique cv segment numbers
ncvseg = length(unique_group_num); 

% Storage for regression vectors
Regcoeff.d = nan(ncvseg,n,ncomp); 
Regcoeff.v = X.v;
Regcoeff.i = 'Regression coefficients per CV segment, Reg(ncvseg, nsamples, ncomp)';


% Cross validation based on id tags (samples with intentical id tags are kept out together)
for k = unique_group_num'
    % Sample indices going into calibration set and out for test set
    indout = (cv_seg ==k); 
    indin = ~indout;

    % PLS calibration based on kept-in samples ----------------------------

    [b0,B,T,W,P,q] = pls_nipals(X.d(indin, :),Y.d(indin,target),ncomp); 

 
    % PLS prediction on kept-out samples ----------------------------------

    YpredTest.d(indout,1:size(B,2)) = b0 + X.d(indout,:)*B;
    YpredTest.i(indout,:) = X.i(indout,:);

    % Regression coefficients ---------------------------------------------
    Regcoeff.d(k,:,1:size(B,2)) = B;
end


% Calculate performance, calibration model on calibration samples ---------
% Build model using all samples.
[b0,B,T,W,P,q] = pls_nipals(X.d,Y.d(:,target),ncomp);
YpredCal.d = b0 + X.d*B;

% Performance for calibration set 
RMSEC = cal_rmse(Y.d(:,target), YpredCal.d);  
R2C = cal_r2(Y.d(:,target), YpredCal.d);

% Calculate performance, calibration model on separate test samples -------
RMSECV = cal_rmse(Y.d(:,target), YpredTest.d);          
RMSECV_corr = rmse_corr(Y.d(:,target), YpredTest.d);
ypred_zerocomp = mean(Y.d(:,target));
rmse_zerocomp = cal_rmse(repelem(ypred_zerocomp,length(Y.d(:,target)))', Y.d(:,target));
R2CV = cal_r2(Y.d(:,target), YpredTest.d);
R2CV_corr = R2_corr(Y.d(:,target), YpredTest.d);
BIAS = cal_bias(Y.d(:,target), YpredTest.d);
slope_intercept = cal_slope(Y.d(:,target), YpredTest.d) ;

% Find optimal num of comps -----------------------------------------------

% Find optimal number of components from given strategy
if  strcmp(p.Results.aoptAlg, 'WM')
    aopt = find_aopt(eval(aoptmetric)); % Based RMSEP/RMSEP_corr, Westad&Martens         
elseif strcmp(p.Results.aoptAlg, 'ALT2')
    aopt = find_aopt([rmse_zerocomp RMSECV], 0.02, 2); % Optimal number of components, Unscr                     
elseif strcmp(p.Results.aoptAlg, 'STATIC')
    aopt = p.Results.ncomp;
end


% Store all important metrics and information -----------------------------

Predictions.residuals = Y.d(:,target) - YpredTest.d;
Predictions.cvpred = YpredTest.d;
Predictions.i = YpredTest.i;
Predictions.v = char(join([repelem(cellstr('#PLS comps: '), ncomp)',cellstr(num2str([1:ncomp]'))]));

Performance.rmsecv = RMSECV; 
Performance.rmsecvcorr = RMSECV_corr;
Performance.rmsec = RMSEC; 
Performance.r2cv = R2CV; 
Performance.r2cvcorr = R2CV_corr; 
Performance.r2c = R2C; 
Performance.biascv = BIAS;
Performance.slopecv = slope_intercept(1); 
Performance.lv = aopt;
Performance.i = ['CV (Leave-idpos-out) ','idpos: ', num2str(idpos)];
Performance.v = char(join([repelem(cellstr('#PLS comps: '), ncomp)',cellstr(num2str([1:ncomp]'))]));     
Performance.regcoeff = B; % Model rebuilt on all samples


% Plots results -----------------------------------------------------------

if strcmp(showPlot,'yes') 
    scrsz = get(0,'ScreenSize');
    
    % RMSE as a function of number of comps, marking the OPTIMAL comp number
    figure('Position',[50 100 scrsz(3)/3 scrsz(4)/2.2])  ; 
    hold on
    pcv = plot(1:ncomp, Performance.rmsecv, '-o','color', plot_colors(1,:));
    pc = plot(1:ncomp, Performance.rmsec, '-o','color', plot_colors(2,:)); 
    box off
    ylabel('RMSE','FontSize', 18)
    xlabel('Num PLS components','FontSize', 18)
    xline(aopt)
    yline(rmse_zerocomp,'--')
    legend('CV', 'Cal')

    % Predicted versus target plot for the OPTIMAL comp number ------------
    figure('Position',[565 100 scrsz(3)/3 scrsz(4)/2.2])  ; 
    hold on
    plot(Y.d(:,target), YpredTest.d(:,aopt),'o', 'color',plot_colors(1,:));    
    if ~isempty(plotlim); xlim(plotlim);ylim(plotlim);end
    text(0.03,0.83,['LV = ', num2str(aopt)],'units','normalized','FontSize',14)
    text(0.03,0.9,['RMSECV = ', num2str(round(Performance.rmsecv(:,aopt),2))],'units','normalized','FontSize',14)
    plot(0:ceil(max(YpredTest.d(:,aopt))), floor(0:ceil(max(YpredTest.d(:,aopt)))));
    xlabel('Target','FontSize',18)
    ylabel('Predicted','FontSize',18)
    box off


end



end

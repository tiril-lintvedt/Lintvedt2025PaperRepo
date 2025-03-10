function X_norm = saisir_normpeak_derivative(X, varargin)
%-------- Normalize spectra by peak ---------------------------------------
%
%   Normalize each spectrum (row) by their value in the x coordinate peakwn
%   Function does not mean center or remove baselines and it may be neccessary  
%   to do this prior to calling this function. 
%
%       INPUT:
%                   X           -      saisir data structure
%                   peakwn      -       x axis value (wavenumber/Raman
%                                       shift) to normalize by.
%
% -------------------------------------------------------------------------
peak_window = 15; % width of peak +- 5 pixels

[nrow ncol]=size(X.d);


if nargin < 3 % Normalise based on one peak -------------------------------
    [~,ipeak] = min(abs(str2num(X.v)- varargin{1}));
    ipeak_lower = ipeak - peak_window;
    ipeak_upper = ipeak + peak_window;
        
    X_norm.d=zeros(nrow,ncol);
    for row= 1:nrow
       if(mod(row,10000)==0)
           disp(num2str([row nrow]));
       end
       
       X_irow = selectrow(X, row); 
       X_sub_deriv1 = saisir_derivative(X_irow,2,9,1); % 1st derivative spectrum
       %saisir_idplot(X_sub_deriv1)
       X_sub_deriv1_peak = selectcol(X_sub_deriv1, ipeak_lower:ipeak_upper);
       peakval = max(X_sub_deriv1_peak.d) - min(X_sub_deriv1_peak.d);       
       X_norm.d(row,:)=X.d(row,:)/peakval; % "OPUS approach"
       
       %X_peak = selectcol(X_irow, ipeak_lower:ipeak_upper);       
       %[X_icorr,baseline,~] = saisir_als(X_peak, 8, 0.01); % correct baseline
       %figure;saisir_idplot(X_icorr); hold on
%        plot(str2num(X_icorr.v), baseline)
%        plot(str2num(X_icorr.v), X_peak.d)
       %[~,ipeak] = min(abs(str2num(X_icorr.v)- varargin{1}));
       %peakval = X_icorr.d(ipeak);
       %X_norm.d(row,:)=X.d(row,:)/peakval;

        
        

    end
    X_norm.v=X.v;
    X_norm.i=X.i;
    
else
    % Normalize based on peak intensity sum -------------------------------
    
    peaks = [];
    for i = 1:(nargin-1) 
        [~,ipeak] = min(abs(str2num(X.v)- varargin{i}));
        peaks(i) = ipeak;
    end   
    
    X_norm.d=zeros(nrow,ncol);   
    
    for row= 1:nrow
       if(mod(row,10000)==0)
           disp(num2str([row nrow]));
       end
       X_norm.d(row,:)=X.d(row,:)/sum(X.d(row,peaks));

    end
    X_norm.v=X.v;
    X_norm.i=X.i;
    
    
end


end
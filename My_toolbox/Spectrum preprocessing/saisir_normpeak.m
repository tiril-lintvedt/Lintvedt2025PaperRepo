function [X_norm, peakval] = saisir_normpeak(X, varargin)
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

[nrow ncol]=size(X.d);
peakval = [];

if nargin < 3 % Normalise based on one peak -------------------------------
    [~,ipeak] = min(abs(str2num(X.v)- varargin{1}));

    X_norm.d=zeros(nrow,ncol);
    for row= 1:nrow
       if(mod(row,10000)==0)
           disp(num2str([row nrow]));
       end
       peakval(row,:) = X.d(row,ipeak);
       X_norm.d(row,:)=X.d(row,:)/X.d(row,ipeak);

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
       peakval(row,:) = sum(X.d(row,peaks));
       X_norm.d(row,:)=X.d(row,:)/sum(X.d(row,peaks));

    end
    X_norm.v=X.v;
    X_norm.i=X.i;
    
    
end


end
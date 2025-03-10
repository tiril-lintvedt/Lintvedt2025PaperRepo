function [Xcorr, baseline] = modpoly(X, varargin)
% [Xcorr, baseline] = modpoly(X, options)
%
% The function modpoly fits and subtracts a polynomial baseline from the
% spectra X. A separate baseline is fitted to each spectrum.
%
% Inputs:
%   - X: A matrix of spectra, dimensions n x m: n samples, m wavelengths
%   - options: A optional struct with the following fields:
%              'order'    : Order of the polynomial.
%                           Default value: 4
%              'interval' : Indices of variables to be used in the fit.
%                           Default: All wavelengths (1:m)
%              'maxIter'  : Maximum iterations per sample.
%                           Default: 1e3.
%              'tolFun'   : Convergence criterium.
%                           Default: 1e-3
%              'plots'    : Turn on (1) or off (0) the final plotting.
%                           Default: 1
% 
% Outputs:
%   - The corrected spectra with baselines subtracted
%   - The baselines that were subtracted
% 
% The algorithm is based on
% Lieber & Mahadevan-Jansen, appl spectrosc. 57 (11), p. 1363 (2003).  
%
% Matlab code written by Nils Kristian Afseth and Martin Høy,
% Nofima Mat AS, Ås, Norway.
% Email: nils.kristian.afseth@nofima.no
%        martin.hoy@nofima.no
%
% February 2009, version 1.0

[nSamp,nVars] = size(X);
meanVal = mean(mean(X));

% Issue warning if X may be transposed
if nSamp > nVars
    warning('modpoly:warn', 'More samples than variables. Should X be transposed?')
end

% Handle required and optional arguments
% First some helper functions to check input arguments
isint    = @(x) isnumeric(x) && isscalar(x) && x>0 && mod(x,1) == 0;
isarray  = @(x) isnumeric(x) && isvector(x) && any(mod(x,1)) == 0 && min(x) > 0 && max(x) <= nVars;
isnumber = @(x) isnumeric(x) && isscalar(x) && x>0;
islogic  = @(x) x == 1 || x == 0;

% Parse input arguments
p = inputParser;
p.addRequired('X',                   @isnumeric);
p.addOptional('order',      4,       isint);
p.addParamValue('interval', 1:nVars, isarray);
p.addParamValue('maxIter',  1e3,     isint);
p.addParamValue('tolFun',   1e-3,    isnumber);
p.addParamValue('plots',    1,       islogic);
p.parse(X, varargin{:});
opts = p.Results;

% Initialise the fitted baseline
baseline = zeros(size(X));
baseline(:,opts.interval) = X(:,opts.interval);

% Loop over all samples
for i=1:nSamp
    for j=1:opts.maxIter

        % Fit a polynomial to each spectrum
        [p, s, mu] = polyfit(opts.interval, baseline(i,opts.interval), opts.order);
        poly = polyval(p, opts.interval, [], mu);

        % The lowest, either the spectrum or the polynomial is selected
        baseline(i,opts.interval) = min( [ baseline(i,opts.interval) ; poly ] );

        % Check for convergence
        if max(abs(baseline(i,opts.interval) - poly)) < opts.tolFun*meanVal
           % disp(sprintf('Sample %2d converged after %3d iterations', i, j));
            break
        end
        
    end
end

% Subtract the fitted baseline from the spectra
Xcorr = X - baseline;

% Optional plotting
if opts.plots
   scrsz = get(0,'ScreenSize');
   figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
   subplot(311)
   plot(X')
   title('Original spectra'), xlabel('Variable index'), ylabel('Signal')
   subplot(312)
   plot(baseline')
   title('Fitted baseline'), xlabel('Variable index'), ylabel('Signal')
   subplot(313)
   plot(Xcorr')
   title('Corrected spectra'), xlabel('Variable index'), ylabel('Signal')
end

function [Y_corr,baseline,wgts] = saisir_als(Ysaisir,lambda,p)
% [Y_corr,baseline,wgts] = saisir_als(Ysaisir,lambda,p)
% -------------------------------------------------------------------------
% Assymetric least squares estimation of baseline, or smoothing if p = 0.5 (no asymmetry) 
% The function uses the SAISIR data structure. For more information:
% http://www.chimiometrie.fr/crbst_15.html#anchor-top 
% -------------------------------------------------------------------------
%   INPUT: 
%           Ysaisir      -      Saisir data structure, spectra
%           lambda       -      scalar, smoothing parameter     
%           p            -      scalar, asymmetric weighting parameter
%
%   OUTPUT: 
%           Y_corr       -      Saisir data structure, the corrected spectra 
%                               (baseline removed / spectra smoothed)  
%           baseline     -      Matrix, baselines for each spectrum
%           wgts         -      Matrix, weights for each spectrum
%
% Reference: Paul Eilers & Hans F. M. Boelens (2005), "Baseline Correction with
% Assymetric Least Squares smoothing)"
%
% Comment: 
% Typical baseline correction parameter values: lambda = 5.8, p = 0.001
% Typical smoothing parameter values: lambda = 0.1,  p = 0.5 (no asymmetry)
%
% -------------------------------------------------------------------------

Y = Ysaisir.d;

% Fit baseline
[n,p1] = size(Y);
baseline = zeros(n,p1);
wgts = zeros(n,p1);
D = diff(speye(p1), 2);
DD = 10^lambda * (D'*D);
for i=1:n
    w = ones(p1,1);
    y = Y(i,:)';
    for it = 1:20 % A fixed maximum number of iterations
        W = spdiags(w, 0, p1, p1);
        C = chol(W + DD);
        z = C \ (C' \ (w .* y));
        w_old = w;
        w = p * (y > z) + (1 - p) * (y < z);
        sw = sum(w_old ~= w);
        if sw == 0, break, end  % Convergence check 
    end
    baseline(i,:) = z;
    wgts(i,:) = w;
end

if p ~= 0.5       % If fitting is asymmetric, do baseline correction
    Y_c = Y - baseline; % Subtract baseline
    Y_c_title = 'Fitted baseline';
    
    
else              % otherwise return smoothed spectrum
    
    Y_c = baseline;
    Y_c_title = 'Smoothed spectrum';
    
end

% Corrected spectra
Y_corr.d = Y_c;
Y_corr.v = Ysaisir.v;
Y_corr.i = Ysaisir.i;


% Control plot
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
subplot(3,1,1)
plot(1:size(Y,2),Y)
ylabel('Original')
xlabel('Variable index')
title('ALS correction','FontSize',18)

subplot(3,1,2)
plot(1:size(Y,2),baseline)
ylabel(Y_c_title)
xlabel('Variable index')

subplot(3,1,3)
plot(1:size(Y,2),Y_corr.d)
ylabel('Corrected spectrum')
xlabel('Variable index')

set(gcf,'Color', [1 1 1])

end
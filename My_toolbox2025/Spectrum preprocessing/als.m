function [Y_c,baseline,wgts] = als(Y,lambda,p)
%  ------------ Asymmetric Least Squares baseline correction --------------
% Assymetric least squares estimation of baseline
%
%   INPUT: 
%           Y            -      Spectra. A column vector, or  matrix with 
%                               nchannels*nspectra matrix
%           lambda       -      scalar, smoothing parameter     
%           p            -      scalar, asymmetric weighting parameter
%
%   OUTPUT: 
%           Y_c          -      Saisir data structure, the corrected spectra 
%                               (baseline removed / spectra smoothed)  
%           baseline     -      Matrix, baselines for each spectrum
%           wgts         -      Matrix, weights for each spectrum
% 
% Ref: Paul Eilers & Hans F. M. Boelens (2005), "Baseline Correction with
% Assymetric Least squares smoothing)"
%
% [baseline,wgts] = als_baseline(y,lambda,p)
%
% Typical baseline correction parameter values: lambda = 5.8, p = 0.001
% Typical smoothing parameter values: lambda = 0.1,  p = 0.5 (no asymmetri)
% -------------------------------------------------------------------------

% Fit baseline
[n,p1] = size(Y);
baseline = zeros(n,p1);
wgts = zeros(n,p1);
D = diff(speye(p1), 2);
DD = 10^lambda * (D'*D);
for i=1:n % for each spectrum fit a curve
    w = ones(p1,1);
    y = Y(i,:)';
    for it = 1:20
        W = spdiags(w, 0, p1, p1);
        C = chol(W + DD);
        z = C \ (C' \ (w .* y));
        w_old = w;
        w = p * (y > z) + (1 - p) * (y < z);
        sw = sum(w_old ~= w);
        if sw == 0, break, end
    end
    baseline(i,:) = z;
    wgts(i,:) = w;
end

if p ~= 0.5       % If fitting is asymmetric, do baseline correction
    % Subtract baseline
    Y_c = Y - baseline;  
    
else              % otherwise return smoothed spectrum
    
    Y_c = baseline;
    
end



% Control plot
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
subplot(3,1,1)
plot(1:size(Y,2),Y)
ylabel('Original')
xlabel('Variable index')

subplot(3,1,2)
plot(1:size(Y,2),baseline)
ylabel('Fitted baseline')
xlabel('Variable index')

subplot(3,1,3)
plot(1:size(Y,2),Y_c)
yline(0);
ylabel('Corrected spectrum')
xlabel('Variable index')

set(gcf,'Color', [1 1 1])

end
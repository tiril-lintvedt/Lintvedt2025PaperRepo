function [SNR_avg, SNR_spec] = cal_SNR(ZSaisir, SG_pol_order, SG_window_size, edge, region)

% ------------------ Signal to Noise calculation --------------------------

% Calculates the SNR for each spectrum in a data set ZSaisir and reports also
% the average SNR in the data set. The methodology is based on
%   a) caluclating a strongly smoothed version of the spectra
%   b) estimating the noise by taking the differenecne between original and
%    strongly smoothed spectrum
%   c) Calculating the standard deviation of the estimated noise.
%   d) Signal is defined as the mean of each spectrum

% The method is similar to the one used in Guo et al. (2020)


% Input:
%           ZSaisir         - Saisir struct, Spectra to calculate SNR for
%           SG_pol_order    - Polynomial order in Savitsky-Golay (SG) smoothing
%           SG_window_size  - Window size in Savitsky-Golay (SG) smoothing
%           edge            - How many point to cut of edges (to avoid SG edge effects)
%           region          - Vector of indices, indicating the wavenumber
%                             regions to use for noise estimation


% Output:
%           SNR_avg         - Number, The average Signal to Noise ratio in data set
%           SNR_spec        - Vector, The SNR for each separate spectrum

% -------------------------------------------------------------------------

if nargin < 2
    SG_pol_order = 2 ;
    SG_window_size = 9;
    edge = 70 ;
    region = 1:length(str2num(ZSaisir.v)); % Full region
end

ZSaisirSmooth = saisir_derivative(ZSaisir,SG_pol_order,SG_window_size,0); % Strongly smoothed spectrum 

% Signal
AverageIntensity = mean(ZSaisir.d,2);

% Estimated absolute noise (standard deviation)
Noise = ZSaisir; 
Noise.d  = ZSaisir.d - ZSaisirSmooth.d; % Isolate noise
Noise = selectcol(Noise, region);  % Use only given region
Noise = selectcol(Noise, edge: length(Noise.v)- edge); % shave off ends to avoid SG edge effects (may not be needed if defined region does not include ends)
% Noise.d = Noise.d(:,edge:(end-edge)); 
% Noise.i = ZSaisir.i;
% Noise.v = ZSaisir.v(edge:(end-edge),:);
Noise_val = std(Noise.d,0,2); % Compute standard deviation along row axis

% Signal to Noise ratio
SNR_spec = AverageIntensity./Noise_val;
SNR_avg =  mean(SNR_spec);

% Control plot, to check that noise isolation and end-shave was successful:
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
subplot(3,1,1)
plot(str2num(ZSaisir.v),ZSaisir.d)
xlim([str2num(ZSaisir.v(1,:)) str2num(ZSaisir.v(end,:))])
title('Noise estimation control check- Absolute noise','FontSize',16)
ylabel('Original', 'FontSize',14)
subplot(3,1,2)
plot(str2num(ZSaisirSmooth.v),ZSaisirSmooth.d)
xlim([str2num(ZSaisir.v(1,:)) str2num(ZSaisir.v(end,:))])
ylabel('Strongly smoothed', 'FontSize',14)
subplot(3,1,3)
plot(str2num(Noise.v), Noise.d)
xlim([str2num(Noise.v(1,:)) str2num(Noise.v(end,:))])
ylabel('Isolated noise','FontSize',14)
xlabel('Raman shift (cm{^{-1}})', 'FontSize',14)
set(gcf, 'Color', [1 1 1])

end
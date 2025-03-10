function Noise = cal_SNR_opus(X, WN1, WN2)
% The OPUS Noise estimation method (for one particular IR application)
%
% -------------------------------------------------------------------------

X1Deriv = saisir_derivative(X,2,61,1); % 1st derivative spectrum % Might 
% need update for derivative hyperparameters , window size , polyn. orddr.

if nargin < 2
    WN1=2000.0; % (choose a region where no chemical bands are expected)
    WN2=2100.0;
end

[y,i1]=min(abs(str2num(X1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(X1Deriv.v)-WN2));
XSaisirWN1_WN2=selectcol(X1Deriv,[i2:i1]); % Might need to change back to i1:i2 
[MaxNoise]=max(XSaisirWN1_WN2.d,[],2);
[MinNoise]=min(XSaisirWN1_WN2.d,[],2);
Noise=MaxNoise-MinNoise;

% Control plot, to check that noise isolation and end-shave was successful
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
subplot(2,1,1)

plot(str2num(X.v), X.d)
title('Noise estimation control check','FontSize',16)
ylabel('Original', 'FontSize',14)
xline(WN1);xline(WN2)

subplot(2,1,2)
plot(str2num(XSaisirWN1_WN2.v), XSaisirWN1_WN2.d)
ylabel('1st Derivative', 'FontSize',14)
xline(WN1);xline(WN2)
set(gcf, 'Color', [1 1 1])


end
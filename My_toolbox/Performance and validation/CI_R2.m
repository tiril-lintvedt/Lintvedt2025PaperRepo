function CI = CI_R2(R2, N)
% Calculation of the double sided 95% confidence interval of the coefficient 
% of determination (R squared) --------------------------------------------
%
%   Input 
%               R2  - Coefficient of determination for a reg. model
%               N   - Number of samples
%
%
% Requirements:
%   * Statistics and Machine learning toolbox
%
% Reference:
% 
% Assumptions: 
%
% -------------------------------------------------------------------------

if nargin == 2
    conf=0.95;  % Default confidence level
end

%N=numel(prediction);
U=conf/2+0.5;
Zeta=0.5*log((1+sqrt(R2))/(1-sqrt(R2))); % Fisher's transformation
CI=[((2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))+1)).^2 , ...
    ((2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))+1)).^2];


% Zeta=0.5*log((1+sqrt(R2))/(1-sqrt(R2))); % Fisher's tranformation
% U = 0.975; % 95% confidence interval
% CI=[((2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))+1)).^2 ,...
%     ((2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))+1)).^2];

end
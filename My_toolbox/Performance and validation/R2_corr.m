function COD = R2_corr(reference,prediction,conf)
% -------- Slope and offset corrected Coefficient of Determination --------

% If prediction has more than one column, separate R2 values are calculated
% for each column of predictions, using the same reference vector.
% Reference must be a one-column vector of reference values.

if nargin == 2
    conf=0.95;
end

%reference=reference(:); prediction=prediction(:);

COD = []; % R2
for i = 1:size(prediction, 2)
    [~,Ri] = MLR(prediction(:,i),reference);
    SSR=sum(Ri.^2);
    SST=sum((reference-mean(reference)).^2);
    COD(i)=1-SSR/SST;
    
end

% Confidence interval
% N=numel(prediction);
% U=conf/2+0.5;
% Zeta=0.5*log((1+sqrt(COD))/(1-sqrt(COD)));
% CI=[((2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))+1)).^2 , ((2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))+1)).^2];
end
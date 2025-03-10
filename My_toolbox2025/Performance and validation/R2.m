function [COD,CI]=R2(reference,prediction,conf)
if nargin == 2
    conf=0.95;
end
reference=reference(:); prediction=prediction(:);
[~,R]=MLR(prediction,reference);

SSR=sum(R.^2);
SST=sum((reference-mean(reference)).^2);
COD=1-SSR/SST;

N=numel(prediction);
U=conf/2+0.5;
Zeta=0.5*log((1+sqrt(COD))/(1-sqrt(COD)));
CI=[((2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta-2*tinv(U,N-1)/sqrt(N))+1)).^2 , ((2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))-1)/(2.718^(2*Zeta+2*tinv(U,N-1)/sqrt(N))+1)).^2];

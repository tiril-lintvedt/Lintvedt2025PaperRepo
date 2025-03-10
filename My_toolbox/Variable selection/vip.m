function VIP = my_vip(W, T, Q, ncomp)
% Variable importance in projections
p = size(W,1);
Q2 = Q(:)'.^2; T2 = sum(T.^2);  % Works for situations where Q can be both a column vector and a row vector
QT = Q2(1:ncomp).*T2(1:ncomp);
WW = W.^2; WW = WW(:,1:ncomp);
WW = WW./sum(WW);
VIP = sqrt(p * sum(WW.*QT,2)./sum(QT,2));   % I changed: sum(QT) - sum(QT,2) for clarity 

% Critical threshold: VIP > 1  (generally used, but noe statistically justified)
end

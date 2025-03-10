function [y,i1] = saisir_getxindex(X, WN)
[y,i1] = min(abs(str2num(X.v)-WN));
end
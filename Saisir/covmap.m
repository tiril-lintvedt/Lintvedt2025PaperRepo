function [saisir] = covmap(saisir1,saisir2)
%covmap 			- assesses the covariances between two tables
%function [cov] = covmap(X1,X2)
%the tables must have the same number of rows
%The results is a p1 x p2 matrix (with p1 and p2 the numbers of rows of X1 and X2 respectively)
%An element cov.d(a,b) is the covariance between the column a of X1 and b of X2


[nrow1 ncol1]=size(saisir1.d);
[nrow2 ncol2]=size(saisir2.d);
if(nrow1~=nrow2); error('The number of rows must be equal');end
a=saisir1.d -ones(nrow1,1)*mean(saisir1.d);
b=saisir2.d -ones(nrow1,1)*mean(saisir2.d);
saisir.d=(1/nrow1)*(a'*b)
saisir.i=saisir1.v; saisir.v=saisir2.v;


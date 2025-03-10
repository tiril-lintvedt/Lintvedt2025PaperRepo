function [saisir] = selectwn(saisir1,wn_region)
%selectcol 		- creates a new data matrix with the selected columns
% usage: [saisir]= selectcol(saisir1,index) 
% saisir is a structure i,v,d 
% the resulting file correspond to the selected columns
% index is a vector of indices (integer) or of booleans

WN1 = wn_region(1);
WN2 = wn_region(end);

[~,i1] = min(abs(str2num(saisir1.v)-WN1));
[~,i2] = min(abs(str2num(saisir1.v)-WN2));

saisir.i=saisir1.i;
saisir.v=saisir1.v(i1:i2,:);
saisir.d=saisir1.d(:,i1:i2);
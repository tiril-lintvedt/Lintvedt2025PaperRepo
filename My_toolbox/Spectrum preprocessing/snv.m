function [X_prep] = snv(X)
% snv 			- Standard normal variate correction on spectra
% function [X1] = snv(X)
% SNV is commonly used in spectroscopy. It basically consists in centering
% and standardizing the ROWS (not the columns) of the data matrix. 
% This procedure may reduce the scatter deformation of spectra
%
% Source: The biospec modelling group, NMBU
%
% Edit : std function ignores nans and inf values


[nrow ncol] = size(X);
X_prep = zeros(nrow,ncol);
X(isinf(X))= nan; % ignore inf values

for row= 1:nrow
   if(mod(row,10000)==0)
       disp(num2str([row nrow]));
   end
   
   
   std1=std(X(row,:), 'omitnan');
   ave=mean(X(row,:),'omitnan');
   X_prep(row,:)=(X(row,:)-ones(1,ncol)*ave)/std1;
end

end

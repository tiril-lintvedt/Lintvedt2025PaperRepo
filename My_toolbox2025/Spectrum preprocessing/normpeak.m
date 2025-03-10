function X_norm = normpeak(wn,X, peakwn)
%-------- Normalize spectra by peak ---------------------------------------
%
%   Normalize each spectrum (row) by their value in the x coordinate peakwn
%   Function does not mean center or remove baselines and it may be neccessary  
%   to do this prior to calling this function. 
%
%       INPUT:
%                   wn          -      number (not char), x axis values
%                   X           -      data array
%                   peakwn      -       x axis value (wavenumber/Raman
%                                       shift) to normalize by.
%
% -------------------------------------------------------------------------

[nrow ncol]=size(X);
[~,ipeak] = min(abs(wn - peakwn));

X_norm=zeros(nrow,ncol);
for row= 1:nrow
   if(mod(row,10000)==0)
       disp(num2str([row nrow]));
   end
   X_norm(row,:)=X(row,:)/X(row,ipeak);
   
end

end
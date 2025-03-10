function mzscore = modified_z_score(x)
% Calculation of the modified Z score for a spectrum x (using median 
% intensity as replacement for mean intensity)
%
%   INPUT:
%               x - a single spectrum
%
%   OUTPUT:
%               mzscore - the modified z score of the spectrum x
%
%   Reference: A Simple Algorithm for Despiking Raman Spectra, Whitaker & Hayes, 2018 
%              python code: https://towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22
%

median_int = median(x); % median intensity in spectrum
mad_int = median(abs(x - median_int)); % MAD = median absolute deviation 
mzscore = 0.6745 * (x - median_int)/mad_int  ; % Modified Z score

end




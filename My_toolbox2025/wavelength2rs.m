function rs = wavelength2rs(wl,laser_wl)
% Calculate Raman shift in wavenumber (cm^-1) from Raman scatter wavelength 

% wl -  Wavelength (nm)
% laser_wl - Laser excitation wavelength (nm)

rs = ((wl/laser_wl)- 1)*((10^7)/wl);

end
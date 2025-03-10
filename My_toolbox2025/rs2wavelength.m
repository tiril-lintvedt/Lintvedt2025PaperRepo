function wl = rs2wavelength(rs,laser_wl)
% Calculate Raman scatter wavelength from Raman shift in wavenumber (cm^-1)  

% rs - Raman shift in wavenumber (cm^-1)
% laser_wl - Laser excitation wavelength (nm)

wl = 1./((1/laser_wl)-(rs/10^7));
end
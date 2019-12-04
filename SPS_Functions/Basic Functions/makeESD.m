% calculates the Energy Spectral Density of a pulse. ESD=abs(fft(y)).^2
% Syntax: [pulse]=makeESD(pulse)
% The function pads the pulse with zeros, shiftes negative frequencies, takes an FFT,
% shifts frequencies back and then crops the data. 
% Added pulse fields are pulse.esd and pulse.esdf

function [freq, amp_spec ] = TakeFFT( t, y )
% This functions performs a one-sided FFT of a function y(t);
%Syntax:
%[frequency, FFT_out] = TakeFFT( t, y )

dt=t(2)-t(1);
Y = fft(y); % Fourier transform of noisy signal
n = size(y,2)/2; % use size for scaling 
amp = abs(Y)/n; % compute amplitude spectrum
f = (0:length(amp)-1)/(2*n*dt); %set frequency
n=floor(length(f)/2); 
freq=f(1:n);         %discard half of the spectrum
amp_spec=amp(1:n);

end


function [pulse]=makeESD(pulse)
% calculates the Energy Spectral Density of a pulse. ESD=abs(fft(y)).^2
% Syntax: [pulse]=makeESD(pulse)
% The function pads the pulse with zeros, shiftes negative frequencies, takes an FFT,
% shifts frequencies back and then crops the data. 
% Added pulse fields are pulse.esd and pulse.esdf
if strcmp(pulse.method, 'Fixed_opBW')==1
    BW=5*pulse.opBW;
else
    BW=100e6;
end
dw=2*BW;

%% pad data to make it a pulse rather than infinite signal
dt=(pulse.t(2)-pulse.t(1));
padn=7.*length(pulse.t);
zeropad=zeros(1, padn);

if size(zeropad,2)~=1
    zeropad=zeropad';
end

t1=0:dt:dt*(padn-1);
t2=(dt*padn)+pulse.t';
t3=t2(end):dt:t2(end)+dt*(padn-1);

tt=[t1, t2, t3];

%% shift negative frequencices to posisitve
fshift=2*abs(min(BW-pulse.f_center));
shifted_signal=pulse.env.*cos(pulse.phi+2*pi*fshift*pulse.t); %shift negative frequencies to positive
padded_sig=[zeropad; shifted_signal; zeropad];

%% take FFT
Fs=1/(pulse.t(2)-pulse.t(1));
L=length(tt);
y=fft(padded_sig);
y=abs(y/L);
Y=y(1:floor(L/2)+1);
Y(2:end-1)=2*Y(2:end-1);
f = Fs*(0:(L/2))/L;
f=f-fshift; %shift back

%% crop output
esdf=f;
esd=Y.^2; %PSD=|f(w)^2|
ind=find(esdf<dw);
esdcrop=esd(ind);
fcrop=esdf(ind);
pulse.esd=esdcrop./max(abs(esdcrop));
pulse.esdf=fcrop;

%check resolution
Fnyq=(fshift+max(pulse.f_mod));
Fs=1/dt;
if  Fs/Fnyq <2
    disp('Warning: Sampling frequency of pulse is not sufficient to accurately caclulate shifted PSD (Fs/F_Nyquist<2). Increase Exp.simpoints')
end

end



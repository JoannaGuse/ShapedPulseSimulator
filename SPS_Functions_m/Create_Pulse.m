
function [pulse] = Create_Pulse(pulse)
%This function generates a fully specified pulse. User must provide B1max,
%Tp and f_max.
%Syntax:
%[ pulse ] = Create_Pulse( pulse )
%pulse is a structure with the following fields:
%
%--------------INPUTS (optional and vary depending on application):-----
% pulse.name: string. pulse name. e.g  'square', 'HSn', 'Hermite180'
% pulse.regime: only for FM pulses. 'Adiabatic', or 'RapidPassage'
% pulse.npts: length of time vector
% pulse.f_center: center frequency (Hz)
% pulse.phi_0: starting phase (radians)
% pulse.Tp: pulse length (seconds) 
% pulse.B1max: B1 strength (Tesla)
% pulse.f_max: Frequency sweep amplitude (Hz)
% Individual pulse parameters e.g. n, beta, sigma
%
%---------------OUTPUTS----------------------
%pulse.class: AM or FM. This dictates how pulse parameters are calculated.
%pulse.Tp: pulse length (seconds) 
%pulse.B1max: B1 strength (Tesla)
%pulse.f_max:  frequency sweep amplitude (Hz)
%pulse.t: pulse time vector (seconds)
%pulse.env: Envelope function i.e. amplitude modulation function (Tesla)
%pulse.f_mod: Frequency modulation function (Hz)
%pulse.phi: Accumulated phase (radians). (i.e. phase modulateion).
%pulse.Bx: pulse.env*cos(pulse.phi) (Tesla)
%pulse.By: pulse.env*sin(pulse.phi) (Tesla)
%pulse.signal: pulse.env*cos(pulse.phi) (Tesla) 
%pulse.esd: Energy Spectral Density of pulse
%pulse.esdf: Frequency axis for ESD (Hz)
%pulse FWHM: FWHM of ESD (linear scale). This is the pulse bandwidth
%pulse.Qcrit: minimum adiabaticity of pulse (only for Adiabatic pulses)
%
%E.g.
%pulse.name='square';
%pulse.B1max=5.6e-4;
%pulse.Tp=32e-9;
%pulse.f_max=0;
%pulse=Create_Pulse(pulse);


%% check inputs
if isfield(pulse, 'name')==0
pulse.name='square';
end

if isfield(pulse, 'method')==0
pulse.method='None';
end

if isfield(pulse, 'B1max')==0
pulse.B1max=5.6e-4;
end

if isfield(pulse, 'Tp')==0
pulse.Tp=32e-9;
end
if isfield(pulse, 'f_max')==0
pulse.f_max=0;
end


if isfield(pulse, 'f_center')==0
    pulse.f_center=0;
end

if isfield(pulse, 'phi0')==0
    pulse.phi0=0;
end

pulse=feval(['Def_', pulse.name], pulse);

pulse.t=linspace(0, pulse.Tp, pulse.npts)';
%% evaluate Bx By and phase
if strcmp(pulse.class, 'AM')==1
    pulse.env=pulse.B1max*pulse.F1;
%   pulse.f_mod=pulse.f_center-0*pulse.t;
 pulse.f_mod=pulse.f_center+0*pulse.t;
    pulse.phi=2*pi*cumtrapz(pulse.t, pulse.f_mod)+pulse.phi0;
elseif strcmp(pulse.class, 'FM')==1
    pulse.env=pulse.B1max*pulse.F1;
%  pulse.f_mod=pulse.f_center-pulse.f_max*pulse.F2;
   pulse.f_mod=pulse.f_center+pulse.f_max*pulse.F2;
    pulse.phi=2*pi*cumtrapz(pulse.t, pulse.f_mod)+pulse.phi0;
else
end

%% make sure outputs are all column vectors;
if size(pulse.env,2)~=1
    pulse.env=pulse.env';
end

if size(pulse.f_mod,2)~=1
    pulse.f_mod=pulse.f_mod';
end
if size(pulse.phi, 2)~=1
    pulse.phi=pulse.phi';
end

pulse.Bx=pulse.env.*cos(pulse.phi);
pulse.By=pulse.env.*sin(pulse.phi);


%% calculate ESD and FWHM
pulse.signal=pulse.env.*exp(1i*pulse.phi);
pulse=makeESD(pulse);

try
pulse.FWHM=FWHM(pulse.esdf,pulse.esd);
end

return



function [pulse]=makeESD(pulse)

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

% %check resolution
% Fnyq=(fshift+max(pulse.f_mod));
% Fs=1/dt;
% if  Fs/Fnyq <2
%     disp('Warning: Sampling frequency of pulse is not sufficient to accurately caclulate shifted PSD (Fs/F_Nyquist<2). Increase Exp.simpoints')
% end

% check to see if there is a sufficient number of points
try
dt=pulse.t(2)-pulse.t(1);
df=1/dt;
middle=round(length(pulse.esd)/2);
cc=pulse.esd(middle:end);
ii=find(abs(cc)<0.1);
fmax=abs(pulse.esdf(middle+ii(1)));
if df/fmax<10
    disp('Insufficent number of timepoints.Errors may occur. Increase pulse.npts.');
end
end

return


function width = FWHM(x,y)
% function width = FWHM(x,y)
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)

% The FWHM result in 'width' will be in units of 'x'

cutoff=0.5; %this is the normalized height at which the width is calculated.

y = y / max(y);
N = length(y);
if y(1) < cutoff                  % find index of center (max or min) of pulse
    [~,centerindex]=max(y);
else
    [~,centerindex]=min(y);
end

%find first crossing
ii = 2;
while sign(y(ii)-cutoff) == sign(y(ii-1)-cutoff)
    ii = ii+1;
end

interp = (cutoff-y(ii-1)) / (y(ii)-y(ii-1));
t1 = x(ii-1) + interp*(x(ii)-x(ii-1));

%start search for next crossing at center
ii = centerindex+1;
while ((sign(y(ii)-cutoff) == sign(y(ii-1)-cutoff)) && (ii <= N-1))
    ii = ii+1;
end

if ii ~= N
    interp = (cutoff-y(ii-1)) / (y(ii)-y(ii-1));
    t2= x(ii-1) + interp*(x(ii)-x(ii-1));
    width = t2- t1;
else
    width = NaN;
end
return

function [F ] = Concatenate_BIR4(pulse)
%Concatenates a pre-defined inversion (pi) pulse into a BIR4 pulse, with rotation angle
%specified by pulse.BIRphase
%e.g.
%pulse.name='HSn'
%pulse.opBW=50e6;
%pulse.theta=pi; 
%pulse=Create_Optimized_Pulse(pulse);
%pulse.BIRphase=pi/3;
%[BIR]=Concatenate_BIR4(pulse);

if isfield(pulse, 'BIRphase')==0
    pulse.BIRphase=pi;
end

n1=round(length(pulse.t)/2);
n2=length(pulse.t)-n1;

AHP.t=pulse.t(1:n1);
AHP.f_mod=pulse.f_mod(1:n1);
AHP.env=pulse.env(1:n1);


RAHP.t=pulse.t(n2:end)-pulse.t(n2);
RAHP.f_mod=pulse.f_mod(n2:end);
RAHP.env=pulse.env(n2:end);


phi1=pi+pulse.BIRphase/2;
phi2=-phi1;

%% make parts
p1=RAHP;
p1.phi=2*pi*cumtrapz(RAHP.t, RAHP.f_mod);

p2=pulse;
p2.phi=2*pi*cumtrapz(p2.t, p2.f_mod)+phi1;

p3=AHP;
p3.phi=2*pi*cumtrapz(AHP.t, AHP.f_mod);

% figure
% plot(p1.t, p1.phi, 'r')
% hold on
% plot(p3.t, p3.phi, 'r')
%% concatenate
F=pulse;

F.env=[p1.env; p2.env; p3.env];
F.f_mod=[p1.f_mod; p2.f_mod; p3.f_mod];
F.phi=[p1.phi; p2.phi+p1.phi(end); p3.phi+p1.phi(end)+p2.phi(end)+phi2];

F.Bx=F.env.*cos(F.phi);
F.By=F.env.*sin(F.phi);
F.t=linspace(0, 2*pulse.t(end), length(F.env))';

F.Tp=F.t(end);

F.signal=F.env.*exp(1i*F.phi);
F=makeESD(F);

try
FWHM = fwhm(F.esdf,F.esd);
F.FWHM=FWHM;
end

end

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

%check resolution
Fnyq=(fshift+max(pulse.f_mod));
Fs=1/dt;
if  Fs/Fnyq <2
    disp('Warning: Sampling frequency of pulse is not sufficient to accurately caclulate shifted PSD (Fs/F_Nyquist<2). Increase Exp.simpoints')
end

end


function [delta,t0] = fwhm(t,y,tr)
%replace NaN with 0
y(isnan(y))=0;


nt = length(y);
dt = diff(t(1:2));
k0 = (1:nt)';
kl = [nt,1:nt-1]';
kr = [2:nt,1]';

if nargin < 3
    [~,ipeak] = max(y);
else
    if length(tr) == 2
        kv = find(tr(1) < t & t < tr(2));
        [~,ipeak] = max(y(kv));
        ipeak = k0(kv(ipeak));
    elseif length(tr) == 1
        itr = interp1(t,k0,tr);
        localmax = find(y(k0) >= y(kl) & y(k0) > y(kr));
        idiff = mod(localmax - itr + nt/2 - 1, nt) - nt/2 + 1;
        [~,i] = min(abs(idiff));
        ipeak = localmax(i);
    end
end

% 3 point quadratic fit to extract peak value

pv = [kl(ipeak),k0(ipeak),kr(ipeak)];
coefs = polyfit((-1:1),[y(pv(1)),y(pv(2)),y(pv(3))],2);
coefs1 = polyder(coefs);

idt = roots(coefs1);
t0 = t(ipeak) + idt*dt;
umax = polyval(coefs,idt);

% find first rising half-max point to left of ipeak

irise = find(y(k0) <= umax/2 & y(kr) > umax/2);
irise = irise + (umax/2 - y(k0(irise)))./(y(kr(irise)) - y(k0(irise)));
irise = min(mod(ipeak-irise,nt));

% find first falling half-max point to right of ipeak

ifall = find(y(k0) >= umax/2 & y(kr) < umax/2);
ifall = ifall + (y(k0(ifall)) - umax/2)./(y(k0(ifall)) - y(kr(ifall)));
ifall = min(mod(ifall-ipeak,nt));

delta = dt*(irise+ifall);
end



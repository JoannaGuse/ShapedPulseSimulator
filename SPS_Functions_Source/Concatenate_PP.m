function [pulse ] = Concatenate_PP(pulse)

% Concatenates a pre-defined pulse into a time-reversed+phase-reversed
% pulse. the output pulse si twice as long as the input pulse
% e.g.
%pulse.name='HSn'
%pulse.opBW=50e6;
%pulse.theta=pi; 
%pulse=Create_Optimized_Pulse(pulse);
%[p2]=Concatenate_PP(pulse);

%% make phase reversal
p1=pulse;

p2=p1;
p2.t=p1.t;
p2.env=flipud(p1.env);
p2.phi=flipud(-p1.phi);
p2.f_mod=flipud(-p1.f_mod);


%% make total
P=pulsetrain1(p1, p2);

pulse.t=P.t;
pulse.Bx=P.Bx;
pulse.By=P.By;
pulse.env=P.env;
pulse.phi=P.phi;
pulse.f_mod=P.f_mod;

[pulse]=makeESD(pulse);
end

function [pulsetrain, train_inputs ] = pulsetrain1( varargin )
%this function makes a pulse train out of input pulses. All info about the
%pulse should be in the pulse input object.
npulses=length(varargin);

%% create individal pulses
train_inputs=cell(npulses);

for jj=1:npulses
    current_pulse=varargin{jj};
    [current_pulse] =Create_Pulse(current_pulse);
    
    if jj==1
        pulsetrain.t=current_pulse.t;
        pulsetrain.Bx=current_pulse.Bx;
        pulsetrain.By=current_pulse.By;
        pulsetrain.env=current_pulse.env;
        pulsetrain.f_mod=current_pulse.f_mod;

    else

        pulsetrain.t=[pulsetrain.t; current_pulse.t(1:end)+pulsetrain.t(end)];
        pulsetrain.Bx=[pulsetrain.Bx; current_pulse.Bx(1:end)];
        pulsetrain.By=[pulsetrain.By; current_pulse.By(1:end)];
        pulsetrain.env=[pulsetrain.env; current_pulse.env(1:end)];
        pulsetrain.f_mod=[pulsetrain.f_mod; current_pulse.f_mod(1:end)];
  
    end
     pulsetrain.phi=2*pi*cumtrapz(pulsetrain.t, pulsetrain.f_mod); 
    train_inputs{jj}=current_pulse;
    clear('current_pulse')
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



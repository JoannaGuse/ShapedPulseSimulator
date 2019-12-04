function pulse = Def_gauss(pulse)
%Gaussian definition function (AM)
%theta: universal
%Shaping parameter: pulse.sigma [0, 1]

pulse.class='AM';
pulse.DynamicVars={'sigma', 'theta'};
pulse.descr='Amplitude modulated pulse for universal rotation operations.';

if isfield(pulse, 'sigma')==0
    pulse.sigma=0.23;
end

if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end
if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

pulse.F1=F1func(pulse);
pulse.BWPfunc=@(pulse)getBWP(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
end

function F1=F1func(pulse)
tau=linspace(0, 1, pulse.npts);
F1=gaussmf(tau, [pulse.sigma,0.5]);
end

function BWP= getBWP(pulse)
BWP=11*exp(-11*pulse.sigma)+1.5*exp(-0.65*pulse.sigma);
BWP=BWP*(-0.0343*pulse.theta^2.315+1.5);
BWP=BWP/1.9;
end


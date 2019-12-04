function pulse = Def_EBURP(pulse)
%EBURP definition function (AM)
%theta:pi/2 (excitation)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for pure-phase Excitation.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end
%add pulse definitions
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 4.84;
pulse.theta=pi/2;

end

function F1=F1func(pulse)
%define F1 (AM function)
tau=linspace(0, 1, pulse.npts);
A0=0.26;
A=[0.91, 0.29, -1.28, -0.05, 0.04, 0.02, 0.06, 0, -0.02];
B=[-0.16, -1.82, 0.18, 0.42, 0.07, 0.07, -0.01, -0.04, 0];
for n=1:length(A)
    term(:,n)=A(n)*cos(2*pi*n*tau)+B(n)*sin(2*pi*n*tau);
end
F1=(A0+sum(term,2));

end
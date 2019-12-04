function pulse = Def_IBURP(pulse)

%IBURP definition function (AM)
%theta:pi (inversion)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for pure-phase Inversion.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 5.21;
pulse.theta=pi;
end

function F1=F1func(pulse)

%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);   
A0=0.5;
A=[0.79, 0, -1.23, -0.19, 0.1, 0.12, 0.04, -0.03, 0.03, -0.01, 0];
B=[-0.71, -1.39, 0.31, 0.47, 0.22, 0.03, -0.05, -0.04, 0, 0.02, 0.01];

for n=1:length(A)
    term(:,n)=A(n)*cos(2*pi*n*tau)+B(n)*sin(2*pi*n*tau);
end
F1=(A0+sum(term,2));

end


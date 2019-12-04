function pulse = Def_REBURP(pulse)
%REBURP definition function (AM)
%theta:pi (universal)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal, pure-phase pi rotation operations.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 5.415;
pulse.theta=pi;

end

function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
A0=0.48;
A=[ -1.03, 1.09, -1.59, 0.86, -0.44, 0.27, -0.17,...
    0.1, -0.08, 0.04, -0.04, 0.01, -0.02, 0, -0.02];
B=zeros(size(A));

for n=1:length(A)
    term(:,n)=A(n)*cos(2*pi*n*tau)+B(n)*sin(2*pi*n*tau);
end
F1=(A0+sum(term,2));

end

function pulse = Def_UBURP(pulse)
%UBURP definition function (AM)
%theta:pi/2 (universal)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal, pure-phase pi/2 rotation operations.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.theta=pi/2;
pulse.BWPfunc=@(pulse) 4.45;

end


function F1=F1func(pulse)

%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
A0=0.27;
A=[ -1.42, -0.33, -1.72, 4.47, -1.33, -0.04,...
    -0.34, 0.5, -0.33, 0.18, -0.21, 0.24, -0.14,...
    0.07, -0.06,+0.06, -0.04, 0.03, -0.03, 0.02];
B=zeros(size(A));

for n=1:length(A)
    term(:,n)=A(n)*cos(2*pi*n*tau)+B(n)*sin(2*pi*n*tau);
end
F1=(A0+sum(term,2));

end
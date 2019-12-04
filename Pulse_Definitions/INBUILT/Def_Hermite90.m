
function pulse = Def_Hermite90(pulse)
%Hermite90 definition function (AM)
%theta:pi/2

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal pi/2 rotation operations.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 3.835;
pulse.theta=pi/2;

end

function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(-1, 1, pulse.npts);
pulse.sigma=0.4;
F1=(1.0-0.667*(tau/pulse.sigma).^2).*exp(-(tau/pulse.sigma).^2);
end
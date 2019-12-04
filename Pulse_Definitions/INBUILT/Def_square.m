function pulse = Def_square(pulse)
%square definition function (AM)
%theta: universal

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal rotation operations.';
pulse.DynamicVars={'theta'};

if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end
if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end


%_______add pulse optimization parameters_________________

pulse.F1=F1func(pulse);
% pulse.BWP=0.8;
pulse.BWPfunc=@(pulse)-2.838e-2*pulse.theta^2.304+1.2;
% pulse.RelArea=1;
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
end


function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
F1=ones(pulse.npts,1);
end







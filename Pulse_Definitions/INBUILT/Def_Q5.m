function pulse = Def_Q5(pulse)
%Q5 definition function (AM)
%theta:pi/2 (universal)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal pi/2 rotation operations.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 6;
pulse.theta=pi/2;

end

function F1=F1func(pulse)

%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
tmax=[16.2, 30.7, 49.7, 52.5, 80.3]./100;
t02=[18.6, 13.9, 14.3, 29.0, 13.7]./2./100;

pulse.A=[-1.48, -4.34, 7.33, -2.30, 5.66];
pulse.sigma=log(2)./(t02).^2;
pulse.tau0=tmax;

for jj=1:5
    g(:,jj)=pulse.A(jj).*exp(-pulse.sigma(jj).*(tau-pulse.tau0(jj)).^2);
    
end
F1=sum(g, 2);
end
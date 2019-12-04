function pulse = Def_G4(pulse)
%G4 definition function (AM)
%theta:pi/2 (excitation)

pulse.class='AM';
pulse.descr='Amplitude modulated Gaussian cascade pulse for pure-phase Excitation.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 5.915;
pulse.theta=pi/2;
end

function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
tmax=[17.7, 49.2, 65.3, 89.2]./100;
t02=[17.2, 12.9, 11.9, 13.9]./100;

pulse.A=[0.62,0.72,-0.91,-0.33];
pulse.sigma=log(2)./(t02).^2;
pulse.tau0=tmax;

for jj=1:4
    g(:,jj)=pulse.A(jj).*exp(-pulse.sigma(jj).*(tau-pulse.tau0(jj)).^2);
    
end
F1=sum(g, 2);

end

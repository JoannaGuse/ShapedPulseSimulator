function pulse = Def_Q3(pulse)
%Q3 definition function (AM)
%theta:pi (universal)

pulse.class='AM';
pulse.descr='Amplitude modulated pulse for universal, pi rotation operations.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 3;
pulse.theta=pi;
end


function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
tmax=[30.6, 54.5, 80.4]./100;
t02=[18.0, 18.3, 24.5]./100;

pulse.A=[-4.39, 4.57, 2.60];
pulse.sigma=log(2)./(t02).^2;
pulse.tau0=tmax;

for jj=1:length(pulse.A)
   g(:,jj)=pulse.A(jj).*exp(-pulse.sigma(jj).*(tau-pulse.tau0(jj)).^2);

end
F1=sum(g, 2);
end

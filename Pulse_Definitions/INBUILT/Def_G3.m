function pulse = Def_G3(pulse)
%G3 definition function (AM)
%theta:pi (inversion)

pulse.class='AM';
pulse.descr='Amplitude modulated Gaussian cascade pulse for pure-phase Inversion.';
pulse.DynamicVars=[];

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

%_______add pulse definitions_________________
pulse.F1=F1func(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);
pulse.BWPfunc=@(pulse) 2.18;
pulse.theta=pi;


end


function F1=F1func(pulse)
%________define F1 (AM function)____________
tau=linspace(0, 1, pulse.npts);
tmax=[28.7, 50.8, 79.5]./100;
t02=[18.9, 18.3, 24.3]./100;

pulse.A=[-1, 1.37, 0.49];
pulse.sigma=log(2)./(t02).^2;
pulse.tau0=tmax;

for jj=1:length(pulse.A)
    g(:,jj)=pulse.A(jj).*exp(-pulse.sigma(jj).*(tau-pulse.tau0(jj)).^2);
    
end
F1=sum(g, 2);

end
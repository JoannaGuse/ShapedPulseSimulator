clear
clc

%define simulation limits
Sys.M0=[0, 0, 1];
Exp.MaxDetuning=80e6;
Exp.npts=201;

%define adiabatic pulse
pulse.name='WURST';
pulse.method='Fixed_B1';
pulse.regime='Adiabatic';
pulse.npts=1e4;
pulse.theta=pi;
pulse.B1max=5e-4; %T
pulse.K=0.9; % bandwidth tuning parameter
pulse.n=10; %shape tuning parameter
pulse= Create_Optimized_Pulse(pulse);

%simulate pulse 
[Sys, Exp, Mx, My, Mz ] = EvolveM_Fast( pulse, Sys, Exp);

%% plot
figure
plot(Exp.detuning/1e6, Mz)
xlabel('Detuning (MHz)')
ylabel('Mz')

 str = {['B_1: ' num2str(pulse.B1max/1e-4, 2),' G'],...
        ['fmax:', num2str(pulse.f_max/1e6, 3)],...
        ['Tp:', num2str(round(pulse.t(end)/1e-9))],...
        ['BW:', num2str(2*pulse.FWHM/1e6, 3)]};
text(0, -0.7, str)
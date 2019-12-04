% 
clear
clc

%define Experimental parameters
Exp.MaxDetuning=40e6;
Exp.npts=200;
Exp.T1=10e-6;
Exp.T2=inf;
Sys.M0=[0,0,1];
Sys.spectrum=ones(1, Exp.npts);
spini=round(Exp.npts/2);

%define pulse sequence
p0.theta=pi/2;
p0.name='square';
p0.method='Fixed_Tp';
p0.Tp=16e-9;
p0.phi0=pi;
p0=Create_Optimized_Pulse( p0 );

wait1=struct('tau', 300e-9, 'npts', 6001);
% simulate
[Sys, Exp, Mx,  My, Mz, Mavg]=SimulatePulseSequence( Sys,Exp, p0, wait1);
t=Exp.t;

Mxavg=Mavg(1,:);
Myavg=Mavg(2,:);
Mzavg=Mavg(3,:);


figure
subplot(1,2,1)
pcolor(t/1e-9, Exp.detuning, My)
shading flat
colorbar
subplot(1,2,2)
plot(t/1e-9, Myavg)
hold on
plot(t/1e-9, My(spini,:))
ylabel('Mz')
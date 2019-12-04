clear
clc
close all

% define experimental parameters
 Exp.MaxDetuning=50e6;
 Exp.npts=301;
 Exp.T1=inf;
 Exp.T2=inf;

% define spin system
Sys.M0=[0,0,1];
Sys.spectrum=ones(1, Exp.npts)';


% define pulse sequence
p1.theta=pi/2;
p1.name='square';
p1=Create_Optimized_Pulse( p1 );

p2.theta=pi;
p2.name='square';
p2=Create_Optimized_Pulse( p2 );

wait1=struct('tau', 600e-9, 'dt', 1e-9); %inter-pulse wait time
wait2=struct('tau', 1200e-9, 'npts', 1200); %detection time

% simulate
[Sys, Exp, Mx,  My, Mz, Mavg]=SimulatePulseSequence( Sys,Exp, p1, wait1, p2, wait2);
t=Exp.t;

Mxavg=Mavg(1,:);
Myavg=Mavg(2,:);
Mzavg=Mavg(3,:);


%plot
figureJo
subplot(2,1,1)
pcolor(t/1e-9,Exp.detuning/1e6, My)
shading flat
title('Hahn Echo')
caxis([-1,1])
colorbar
ylabel('Detuning (MHz)')
subplot(2,1,2)
plot(t/1e-9, Myavg, 'b')
xlabel('Time (ns)')
ylabel('My_{avg}')

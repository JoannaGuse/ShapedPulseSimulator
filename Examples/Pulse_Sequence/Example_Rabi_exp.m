
clear
clc

% define experimental parameters
Exp.MaxDetuning=80e6;
Exp.npts=301;
Exp.T1=100e-6;
Exp.T2=500e-9;

Sys.M0=[0,0,1];
Sys.spectrum=ones(1, Exp.npts)';

% define pulse sequence
p2.theta=pi;
p2.name='square';
p2.Tp=32e-9;
p2.method='Fixed_Tp';
p2=Create_Optimized_Pulse( p2 );

wait1=struct('tau', 250e-9, 'npts', 501);  %inter-pulse delay
wait2=struct('tau', 600e-9, 'npts', 1201); % detection time 

% define nutation experiment
param.d30=2e-9;  %pulse length increment
param.sx=51;     %number of increments (size of x axis)
param.pg=150e-9; % echo 'integration' width
plengths=0:param.d30:param.d30*param.sx; %vector of nutation pulse lengths
tau=plengths/1e-9;

% set up loop
hh=waitbar(0, 'looping over position displacement');
Sys0=Sys;
for jj=1:length(plengths);
clear Sys p1
p1=p2;
p1.Tp=plengths(jj);
p1=Create_Pulse(p1);
Sys=Sys0;

% simulate
[Sys, Exp, Mx,  My, Mz, Mavg]=SimulatePulseSequence( Sys,Exp, p1, wait1, p2, wait2);
t=Exp.t;

MY(:,jj)=Mavg(2,:);
tt(:,jj)=t;

%find top of echo and integrate
d0=p1.t(end)+p2.t(end)+2*wait1.tau; 
width=param.pg;
t0=d0-width/2;
t1=d0+width/2;
echostart(jj)=t0;
echoend(jj)=t1;

[~, i1]=min(abs(t-t0));
[~, i2]=min(abs(t-t1));
echo=MY(i1:i2, jj);
techo=t(i1:i2);
intecho(jj)=sum(echo);
ech(:,jj)=echo';
tech(:,jj)=techo';


waitbar(jj/length(plengths), hh)
end
close(hh)


%% plot

figureJo
clims=[-0.3, 0.3];
subplot(1,3,[1,2])
pcolor(tt',tau', MY')
shading flat
colorbar
title(['My ', p1.name, '\theta - ', num2str(round(p2.Tp/1e-9)),' ns ', p2.name ' \pi'])
ylabel('\tau (ns)')
xlabel('time (ns)')
caxis(clims)

ax1 = gca;
ax2 = axes('position',get(ax1,'position'),'color','none');
set(ax2,'YAxisLocation','right', 'XAxisLocation','top')
hold on
plot( echostart, tau,'color', 'w', 'linewidth', 2)
plot(echoend, tau, 'color', 'w', 'linewidth', 2)
set(ax2,'ylim',[tau(1), tau(end)]);
set(ax2,'Visible','off');
linkaxes([ax1 ax2], 'x');


subplot(1,3,3)
plot(tau, intecho, 'k')
xlabel('\tau (ns)')
ylabel('integrated echo')
hold on





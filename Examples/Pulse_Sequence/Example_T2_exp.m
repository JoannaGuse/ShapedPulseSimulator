
clear
clc
close all

% define experimental parameters
Exp.MaxDetuning=80e6; %choose carefully. The echo is summed over the these detunings. MaxDetuning>min[pulse BW,Sys.spectrum extent]
Exp.npts=301; % make sure you have enough points so echo signal is an accurate average
Exp.T1=100e-6;
Exp.T2=500e-9;

Sys.M0=[0,0,1];
Sys.spectrum=ones(1, Exp.npts);

% define pulse sequence
p1.theta=pi/2;
p1.name='square';
p1=Create_Optimized_Pulse( p1 );
p1.invBW=round(p1.FWHM/1e6)*1e6;

p2.theta=pi;
p2.name='square';
p2=Create_Optimized_Pulse( p2 );

wait1=struct('tau', 250e-9, 'npts', 501);
wait2=struct('tau', 2000e-9, 'npts', 2001);

param.pg=5*1/(p1.invBW); %width of Hahn echo
param.d30=25e-9; %pulse position increment
param.sx=21;     %number of increments


% set up loop
positions=0:param.d30:param.d30*param.sx; %vector of tau increments
tau0=wait1.tau; %initial tau
tau=tau0+positions; 

hh=waitbar(0, 'looping over position displacement');
Sys0=Sys;
for jj=1:length(positions);
wait1.tau=tau0+positions(jj);
clear Sys
Sys=Sys0;
% simulate
[Sys, Exp, Mx,  My, Mz, Mavg]=SimulatePulseSequence( Sys,Exp, p1, wait1, p2, wait2);
t=Exp.t;
MY(:,jj)=Mavg(2,:);
tt(:,jj)=t;

%find top of echo and integrate
d0=p1.t(end)+p2.t(end)+2*wait1.tau; 
width=5*1/(p1.invBW);
t0=d0-width/2;
t1=d0+width/2;

[~, i1]=min(abs(t-t0));
[~, i2]=min(abs(t-t1));
echo=MY(i1:i2, jj);
techo=t(i1:i2);
intecho(jj)=sum(echo);
ech(:,jj)=echo';
tech(:,jj)=techo';

waitbar(jj/length(positions), hh)
end
close(hh)

T=2*positions; %correct tau axis is 2*increment lenght

%% plot
figureJo
subplot(1,3,1)
pcolor(tau/1e-9,tt/1e-9,MY)
shading flat
xlabel('\tau (ns)')
ylabel('time (ns)')
subplot(1,3,2)
plot(tech, ech)
xlabel('time (ns)')
ylabel('My')
hold on
subplot(1,3,3)

plot(T, intecho, 'k')
hold on
plot(T, max(max(abs(intecho)))*exp(-T*(1/Exp.T2+1/Exp.T1)), 'g--')
legend('sim', '|intEcho(0)|e^{-tau"(1/T_2+1/T_1)}') %note: 2*tau gives the time. 
xlabel('\tau (ns)')
ylabel('integrated echo')

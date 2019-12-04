clear
clc
close all

addpath(genpath(pwd))
addpath(genpath('C:\Users\Joanna\Dropbox (UNSW ESR Lab)\Joanna\PULSE DESIGN'))

%% define experimental parameters
 Exp.MaxDetuning=50e6;
 Exp.npts=201;
 Exp.T1=100e-6;
 Exp.T2=800e-9;

% define spin system
Sys.M0=[0,0,1];
Sys.spectrum=ones(1, Exp.npts)';

% Exp.detuning=linspace(-Exp.MaxDetuning, Exp.MaxDetuning, Exp.npts);
% sigma=30e6;
% Sys.spectrum=gaussmf(Exp.detuning, [sigma, 0]);

%% define pulse sequence   INV-t1-pi/2-tau-pi
pINV.theta=pi;
pINV.name='square';
pINV.method='Fixed_Tp';
pINV.Tp=32e-9;
pINV=Create_Optimized_Pulse( pINV );

p2.theta=pi/2;
p2.name='square';
p2.method='Fixed_Tp';
p2.Tp=16e-9;
p2=Create_Optimized_Pulse( p2 );

p3.theta=pi;
p3.name='square';
p3.method='Fixed_Tp';
p3.Tp=32e-9;
p3=Create_Optimized_Pulse( p3 );

wait1=struct('tau', 50e-9, 'npts', 100);
wait2=struct('tau', 250e-9, 'npts', 250);
wait3=struct('tau', 600e-9, 'npts', 600);

pg=200e-9; %echo width 
d0=wait1.tau+2*wait2.tau+pINV.t(end)+p2.t(end)+p3.t(end)-pg/2; %echo start 
d1=d0+pg; %echo end

ifplot=0;

%% simulate
Sys0=Sys;

fsweep=linspace(-50e6, 50e6, 51); %center frequency sweep position of 'ELDOR' pulse
hh=waitbar(0, 'looping over f');

if ifplot==1
figureJo;
end

for jj=1:length(fsweep)  
    clear Sys
    Sys=Sys0;
    
%increment center frequency position of inversion pulse
f_ELDOR=fsweep(jj);
p1=pINV;
p1.f_center=f_ELDOR;
p1=Create_Pulse( p1 );

% simulate
[Sys, Exp, Mx,  My, Mz, Mavg]=SimulatePulseSequence( Sys,Exp, p1, wait1, p2, wait2, p3, wait3);
t=Exp.t;

Mx1(:,jj)=Mavg(1,:);
My1(:,jj)=Mavg(2,:);
Mz1(:,jj)=Mavg(3,:);

[~,m]=min(abs(t-d0));
[~,n]=min(abs(t-d1));


if ifplot==1
subplot(1,2,1)
pcolor(t, Exp.detuning, My)
shading flat
subplot(1,2,2)
plot(t, My1(:,jj), 'b')
hold on
plot(t(m:n), My1(m:n, jj), 'r')
end

waitbar(jj/length(fsweep), hh)
end
close(hh)

% caclulate inversion spectrum
[~,m]=min(abs(t-d0));
[~,n]=min(abs(t-d1));
echo=mean(My1(m:n, :));
maxecho=max(echo);
invspec=echo./maxecho;

%% plot
figure('color', [1,1,1], 'units', 'normalized', 'position', [.1,.1,.7,.4])
fa=16;
subplot(1,3,1)
pcolor(t/1e-9,fsweep/1e9, Mx1')
shading flat
xlabel('Time (ns)', 'Fontsize', fa)
ylabel('Freq (GHz)', 'Fontsize', fa)
caxis([-.2, 0.2])
colorbar
set(gca, 'Fontsize', fa)
title(['Mx ', p1.name, p2.name, p3.name], 'Fontsize', fa)
subplot(1,3,2)
pcolor(t/1e-9,fsweep/1e9, My1')
shading flat
xlabel('Time (ns)', 'Fontsize', fa)
ylabel('Freq (GHz)', 'Fontsize', fa)
caxis([-.2, 0.2])
colorbar
title('My', 'Fontsize', fa)
set(gca, 'Fontsize', fa)

subplot(1,3,3)
plot(fsweep/1e6, invspec,'b', 'linewidth', 2)
axis tight
ylim([-1, 1])
ylabel(' My_{avg}', 'Fontsize', fa)
xlabel('Freq (GHz)', 'Fontsize', fa)
set(gca, 'Fontsize', fa)


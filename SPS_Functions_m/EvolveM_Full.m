function [Sys, Exp, Mx, My, Mz]=EvolveM_Full(pulse, Sys, Exp)
%This function solves the Bloch equations (with T1 and T2 relaxation) to
%simulate the evolution of an ensemble of magnetization vectors subject to 
%a time-dependent pulse. The function returns the full magnetization M(time, detuning).
%
%syntax:[ Sys,Exp, Mx, My, Mz ] = EvolveM_Fast( pulse, Sys, Exp);
%
%________________INPUTS_____________________
%pulse                :pulse object which contains fields t, Bx(t) and By(t) 
%
%Sys.M0               :initial magnetization. Either [1, 3] or [Exp.npts,3] matrix. Defualt: [0,0,1];
%
%Exp.npts             :number of x axis points(detuning values).Defualt:101
%Exp.MaxDetuning      :x axis detuning limit (Hz). Default:100e6
%Exp.T1               :T1 relaxation time (seconds). Default:inf
%Exp.T2               :T2 relaxation time (seconds). Default:inf
%Exp.gamma            :gyromagnetic ratio (%rad/s/T). Default: 1.760859e11
%
%_____________OUTPUTS (and additional fields)____
%Mx, My, Mz           :Magnetization at each detuning value during the pulse.
%Sys.M0               :Full initial magnetization. [Exp.npts, 3] matrix
%Exp.detuning         :x axis detuning frequencies (Hz) defined between
%                     [-Exp.MaxDetuning, Exp.MaxDetuning ]
%Exp.Bz               :Bz component of the effective B field in the rotating
%                      frame.
%
%Example: 
%Sys.M0=[0,0,1];
%Exp.T1=1e6;
%Exp.T2=500e-9;
%Exp.npts=501;
%Exp.MaxDetuning=100e6;
%[ Sys,Exp, Mx, My, Mz ] = EvolveM_Full( pulse, Sys,Exp);(use defualt values)

%% check inputs
if isfield(Exp, 'gamma')==0 || isempty(Exp)==1
   Exp.gamma=1.760859e11; %rad/s/T %DO NOT CHANGE
end

if isfield(Exp, 'npts')==0
    Exp.npts=101;
end

if isfield(Exp, 'MaxDetuning')==0
    Exp.MaxDetuning=100e6;
end

if isfield(Exp, 'T1')==0
    Exp.T1=inf;
end
if isfield(Exp, 'T2')==0
    Exp.T2=inf;
end

if isfield(Exp, 'Mz_eqilib')==0
    Exp.Mz_equilib=1;
end

if exist('Sys', 'var')==0 || isempty(Sys)==1
    Sys.M0=[0,0,1];
end

if size(Sys.M0, 2)==3 && size(Sys.M0, 1)==1
Sys.M0=repmat(Sys.M0, Exp.npts, 1);
elseif size(Sys.M0, 1)==3 && size(Sys.M0, 2)==1
Sys.M0=repmat(Sys.M0', Exp.npts, 1);
elseif size(Sys.M0, 2)==3 && size(Sys.M0, 1)==Exp.npts
Sys.M0=Sys.M0;
elseif size(Sys.M0, 1)==3 && size(Sys.M0, 2)==Exp.npts
Sys.M0=Sys.M0';
else
    disp('incorrect Sys.M0 size')
end


%% perform evolution
dB=Exp.MaxDetuning*2*pi/Exp.gamma;
Exp.Bz=linspace(-dB, dB, Exp.npts);
Exp.detuning=Exp.Bz*Exp.gamma/(2*pi);

%simulate spins
Mtot=EvolveM(pulse, Sys, Exp);
Mx=squeeze(Mtot(:, 1,:)); 
My=squeeze(Mtot(:, 2,:)); 
Mz=squeeze(Mtot(:, 3,:));

return


function  [Mtot]=EvolveM(pulse, Sys, Exp)
dt=pulse.t(2)-pulse.t(1);
M = Sys.M0;  %this sets up an inital magnetization direction for all the spins.

%set up storage matrices
Mtot=zeros(size(M, 1),size(M, 2), length(pulse.t));

%set up B field sizes
B0 = zeros(length(Exp.Bz),3);
B0(:,3)= Exp.Bz;        
B1 = zeros(length(Exp.Bz),3);

Bx=pulse.Bx+1e-15;
By=pulse.By+1e-15; %add offset so we don't divide by 0

for tt=1:length(pulse.t)
    B1(:,1)=Bx(tt);
    B1(:,2)=By(tt);
    Beff = B0 + B1;
  
    M= EvolveSpins2(M,Beff,repmat(dt,length(Exp.detuning),1), Exp); 
    Mtot(:,:,tt)=M;
end


function e = EvolveSpins2(e,Be,tevolve, Exp)

% Exp.gamma=1.7608597e11; %electron gyromagnetic ratio
ne = New_Div(Be,sqrt(sum(Be.^2,2)));
Be = sqrt(sum(Be.^2,2));
e = Rotation(e,ne,-Exp.gamma*New_Mult(Be,tevolve));
e=Relax(e,tevolve(1),Exp);
return

function [M]=Relax(a,dt,Exp)
Mx=a(:,1)*exp(-dt/Exp.T2);
My=a(:,2)*exp(-dt/Exp.T2);
Mz=Exp.Mz_equilib-(Exp.Mz_equilib-a(:,3))*exp(-dt/Exp.T1);

M=[Mx, My, Mz];
return

%%_______________________________________________
function output = Rotation(r,n,theta)
output = New_Mult(r,cos(theta)) + New_Mult(New_Mult(n,(sum(n.*r,2))),(1-cos(theta))) + New_Mult(cross(r,n),sin(theta));
return
%%_______________________________________________
function output = New_Mult(a,b)
output = a;
for n = 1:length(b)
    output(n,:) = output(n,:)*b(n);
end
return

%%_______________________________________________

function output = New_Div(a,b)
output = a;
for n = 1:length(b)
    output(n,:) = output(n,:)/b(n);
end
return


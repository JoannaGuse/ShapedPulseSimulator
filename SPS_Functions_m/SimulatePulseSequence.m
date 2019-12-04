function [Sys, Exp, Mx, My, Mz, Mavg ] = SimulatePulseSequence(Sys, Exp, varargin )
 
%simulates pulse sequence using Bloch equations.
%[Sys, Exp, Mx, My, Mz, Mavg ] = SimulatePulseSequence(Sys, Exp, varargin);
%
%________________INPUTS_____________________
%
%Exp.npts             :number of x axis points(detuning values).Defualt:101
%Exp.MaxDetuning      :x axis detuning limit (Hz). Default:100e6
%Exp.T1               :T1 relaxation time (seconds). Default:inf
%Exp.T2               :T2 relaxation time (seconds). Default:inf
%Exp.gamma            :gyromagnetic ratio (%rad/s/T). Default: 1.760859e11
%
%Sys.M0               :initial magnetization.[Exp.npts, 3] matrix
%Sys.spectrum         :spin spectrum. Used to weight the average magnetization output. 
%
%pulse                :pulse object which contains t, Bx(t) and By(t) 
%
%wait.tau             :free evolution time (seconds)                
%wait.dt              :free evolution time-step (seconds)(optional)
%wait.npts            :number of free evolution timepoints (optional)
%
% 
%_____________OUTPUTS (and additional fields)____
%Mx, My, Mz           :Magnetization matrices at each detuning value and time.
%Mavg                 :Average magnetization [3, t], weighted by the spectrum
%Sys.spectrum         :Spin spectrum. Used to weight the magnetization output.
%Exp.detuning         :x axis detuning frequencies (Hz) defined between
%                     [-Exp.MaxDetuning, Exp.MaxDetuning ]
%Exp.Bz               :Bz component of the effective B field in the rotating frame (Tesla)
%Exp.t                :pulse sequence concatenated time vector (seconds)
%
%Example:
%
%Sys.M0=[0,0,1];
%Exp.T2=500e-9;
%p1=struct('name', 'square','theta',pi/2','opBW', 50e6);
%p1=Create_Optimized_Pulse(p1);
%wait1=struct('tau', 250e-9, 'dt', 1e-9);
%p2=struct('name', 'square','theta',pi','opBW', 50e6);
%p2=Create_Optimized_Pulse(p2);
%wait2=struct('tau', 500e-9, 'dt', 1e-9);
%[Sys, Exp, Mx, My, Mz, Mavg ] = SimulatePulseSequence(Sys, Exp, p1,wait1, p2, wait2);
%plot(Exp.t, Mavg(2,:));

 %% check inputs
if isfield(Exp, 'gamma')==0 || isempty(Exp)==1
   Exp.gamma=1.760859e11; %rad/s/T %DO NOT CHANGE
end
 
if isfield(Exp, 'MaxDetuning')==0
    Exp.MaxDetuning=100e6;
end
 
 
if isfield(Exp, 'npts')==0
    Exp.npts=round(Exp.MaxDetuning/0.25e6);
end

Resolution=Exp.MaxDetuning/Exp.npts;
if Resolution>0.5e6
    disp('insufficient Exp.npts. may result in distorted spectra')
    Exp.npts=round(Exp.MaxDetuning/0.25e6);
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
 
dB=Exp.MaxDetuning*2*pi/Exp.gamma;
Exp.Bz=linspace(-dB, dB, Exp.npts);
Exp.detuning=Exp.Bz*Exp.gamma/(2*pi);

%% start evolution
p=varargin;
Mstart=Sys.M0;
Mx=[];
My=[];
Mz=[];
t=0;

for jj=1:length(varargin)   
    if isfield (p{jj}, 'name')==1
%         disp(p{jj}.name)
        M=EvolveM(p{jj}, Sys, Exp);
        M=squeeze(M);  
        if isvector(Sys.M0)
            Mx=[Mx,M(1,:)];
            My=[My,M(2,:)];
            Mz=[Mz,M(3,:)];
            Mend=M(:,end);
            t=[t; t(end)+p{jj}.t];
            Sys.M0=Mend;
        else
            Mx=[Mx, squeeze(M(:,1,:))];
            My=[My, squeeze(M(:,2,:))];
            Mz=[Mz, squeeze(M(:,3,:))];
            t=[t; t(end)+p{jj}.t];
            Mend=squeeze(M(:,:,end));
            Sys.M0=Mend;
          
        end
    else 
        if isfield(p{jj}, 'npts')==1
          tevol=linspace(0, p{jj}.tau, p{jj}.npts)';
        elseif isfield(p{jj}, 'dt')==1
          tevol=(0:p{jj}.dt:p{jj}.tau)';
        else
          tevol=(0:1e-9:p{jj}.tau)';
        end
        [~, ~, M]=Free_precess(Sys, Exp, tevol);
        M=squeeze(M);     
        if isvector(Sys.M0)
            Mx=[Mx,M(1,:)];
            My=[My,M(2,:)];
            Mz=[Mz,M(3,:)];
            Mend=M(:,end);
            t=[t; t(end)+tevol];
            Sys.M0=Mend';
        else
            Mx=[Mx, squeeze(M(:,1,:))];
            My=[My, squeeze(M(:,2,:))];
            Mz=[Mz, squeeze(M(:,3,:))];
            t=[t; t(end)+tevol];
            Mend=squeeze(M(:,:,end));
            Sys.M0=Mend;
        end
       
    end
end
Exp.t=t;
Mx=[Mstart(:,1), Mx];
My=[Mstart(:,2), My];
Mz=[Mstart(:,3), Mz];

if isfield(Sys, 'spectrum')==0
     Sys.spectrum=ones(size(Exp.Bz))';
end

if numel(Sys.spectrum)~=length(Sys.spectrum)
    disp('Incorrect spectrum size. Default to uniform spectrum')
    Sys.spectrum=ones(size(Exp.Bz))';
    gg=repmat(Sys.spectrum,1, length(t));     
elseif size(Sys.spectrum,1)==1
    Sys.spectrum=Sys.spectrum';
    gg=repmat(Sys.spectrum,1, length(t));
else   
    gg=repmat(Sys.spectrum,1, length(t));
end

Mx=Mx.*gg;
My=My.*gg;    
Mz=Mz.*gg;

Mxavg=mean(Mx, 1)/sum(Sys.spectrum)*Exp.npts;
Myavg=mean(My, 1)/sum(Sys.spectrum)*Exp.npts;
Mzavg=mean(Mz, 1)/sum(Sys.spectrum)*Exp.npts;

Mavg=[Mxavg; Myavg; Mzavg];

return


function [Sys, Exp, M]=Free_precess(Sys, Exp, t)
% Perform Free evolution during time t
% [M]=Free_precess(Sys, Exp, t);
% Sys.M0          :[3, Exp.npts] matrix of initial Magnetizations.
% Exp.T1, Exp.T2  :T1 and T2 rlaxation rates
% Exp.detuning    :[1, Exp.npts] vecotr of detuning values (Hz)
% t               :Free evolution time vector (s) 

dt=t(2)-t(1);
M=zeros(Exp.npts, 3, length(t));
for jj=1:Exp.npts   
    df=Exp.detuning(jj);
    [A,B] = freeprecess1(dt,Exp,df);       
    Mj = zeros(3,length(t));	% Keep track of magnetization at all time points.
    Mj(:,1)=Sys.M0(jj,:)';
    
    for k=2:length(t)
        Mj(:,k) = A*Mj(:,k-1)+B;
    end;
    M(jj,:,:)=Mj; 
  
end

return

function [Afp,Bfp]=freeprecess1(T,Exp,df)
phi = 2*pi*df*T;	% Resonant precession, radians.
E1 = exp(-T/Exp.T1);	
E2 = exp(-T*(1/Exp.T1+1/Exp.T2));
Rz = [cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];

Afp = [E2 0 0;0 E2 0;0 0 E1]*Rz;
Bfp = [0 0 1-E1]';
return

%______________________________________________________
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
ne = New_Div(Be,sqrt(sum(Be.^2,2)));
Be = sqrt(sum(Be.^2,2));
e = Rotation(e,ne,-Exp.gamma*New_Mult(Be,tevolve));
e=Relax(e,tevolve(1),Exp);
return
 
function [M]=Relax(a,dt,Exp)

Mx=a(:,1)*exp(-dt/Exp.T2);
My=a(:,2)*exp(-dt/Exp.T2);
Mz=Exp.Mz_equilib-(Exp.Mz_equilib-a(:,3))*exp(-dt/Exp.T1);

% Mzeq=1;
% T1=inf;
% T2=inf;
% Mx=a(:,1)*exp(-dt/T2);
% My=a(:,2)*exp(-dt/T2);
% Mz=Mzeq-(Mzeq-a(:,3))*exp(-dt/T1);

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
function [Sys, Exp, Mx, My, Mz]=Free_precess(Sys, Exp, t)
% Performs free evolution during time t, including T1 and T2 relaxation.
% [Sys, Exp, Mx, My, Mz]=Free_precess(Sys, Exp, t)
%
% _________Inputs_______
% Sys.M0          :[3, Exp.npts] matrix of initial magnetization vectors.
% Exp.T1          :T1 relaxation rate (seconds)
% Exp.T2          :T2 relaxation rate (seconds)
% Exp.npts        :number of detuning values
% Exp.detuning    :[1, Exp.npts] vector of detuning values (Hz)
% t               :Free evolution time vector (seconds) 
%
%_________Outputs_______
%Mx, My, Mz       :Magnetization at each detuning value during free evolution.
%
%Example
% 
% Exp.T1=1e-6;
% Exp.T2=500e-9;
% Exp.npts=101;
% Exp.detuning=linspace(-50e6, 50e6, Exp.npts);
% t=linspace(0, 500e-9, 501);
% Sys.M0=repmat([1,0,0], Exp.npts, 1);
% [Sys, Exp, Mx, My, Mz]=Free_precess(Sys, Exp, t);
% plot(t, Mx)

Exp.npts
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
  Mx=squeeze(M(:,1,:));
  My=squeeze(M(:,2,:));
  Mz=squeeze(M(:,3,:));
return

function [Afp,Bfp]=freeprecess1(T,Exp,df)
phi = 2*pi*df*T;	% Resonant precession, radians.
E1 = exp(-T/Exp.T1);	
E2 = exp(-T*(1/Exp.T1+1/Exp.T2));
Rz = [cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];

Afp = [E2 0 0;0 E2 0;0 0 E1]*Rz;
Bfp = [0 0 1-E1]';
return





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

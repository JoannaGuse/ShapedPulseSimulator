
%This function solves the Bloch equations (with T1 and T2 relaxation) to
%simulate the evolution of an ensemble of magnetization vectors subject to 
%a time-dependent pulse. The function returns the end-of-pulse magnetization M(Tp, detuning).
%
%syntax:[ Sys,Exp, Mx, My, Mz ] = EvolveM_Full( pulse, Sys, Exp);
%________________INPUTS_____________________
%pulse                :pulse object which contains fields t, Bx(t) and By(t) 
%Sys.M0               :initial magnetization. Either [1, 3] or [Exp.npts,3] matrix. Defualt: [0,0,1];
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
%Examples:  
%[ Sys,Exp, Mx, My, Mz ] = EvolveM_Fast( pulse, [],[]);(use defualt values)
%
%Sys.M0=[0,0,1];
%Exp.T1=1e6;
%Exp.T2=500e-9;
%Exp.npts=501;
%Exp.MaxDetuning=100e6;
%[ Sys,Exp, Mx, My, Mz ] = EvolveM_Fast( pulse, Sys,Exp);(use defualt values)
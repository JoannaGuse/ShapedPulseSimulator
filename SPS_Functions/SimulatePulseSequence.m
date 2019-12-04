 
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
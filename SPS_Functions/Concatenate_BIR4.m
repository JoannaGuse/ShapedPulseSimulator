
%Concatenates a pre-defined inversion (pi) pulse into a BIR4 pulse, with rotation angle
%specified by pulse.BIRphase
%e.g.
%pulse.name='HSn'
%pulse.opBW=50e6;
%pulse.theta=pi; 
%pulse=Create_Optimized_Pulse(pulse);
%pulse.BIRphase=pi/3;
%[BIR]=Concatenate_BIR4(pulse);

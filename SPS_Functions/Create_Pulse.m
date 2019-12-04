%This function generates a fully specified pulse. User must provide B1max,
%Tp and f_max.
%Syntax:
%[ pulse ] = Create_Pulse( pulse )
%pulse is a structure with the following fields:
%
%--------------INPUTS (optional and vary depending on application):-----
% pulse.name: string. pulse name. e.g  'square', 'HSn', 'Hermite180'
% pulse.regime: only for FM pulses. 'Adiabatic', or 'RapidPassage'
% pulse.npts: length of time vector
% pulse.f_center: center frequency (Hz)
% pulse.phi_0: starting phase (radians)
% pulse.Tp: pulse length (seconds) 
% pulse.B1max: B1 strength (Tesla)
% pulse.f_max: Frequency sweep amplitude (Hz)
% Individual pulse parameters e.g. n, beta, sigma
%
%---------------OUTPUTS----------------------
%pulse.class: AM or FM. This dictates how pulse parameters are calculated.
%pulse.Tp: pulse length (seconds) 
%pulse.B1max: B1 strength (Tesla)
%pulse.f_max:  frequency sweep amplitude (Hz)
%pulse.t: pulse time vector (seconds)
%pulse.env: Envelope function i.e. amplitude modulation function (Tesla)
%pulse.f_mod: Frequency modulation function (Hz)
%pulse.phi: Accumulated phase (radians). (i.e. phase modulateion).
%pulse.Bx: pulse.env*cos(pulse.phi) (Tesla)
%pulse.By: pulse.env*sin(pulse.phi) (Tesla)
%pulse.signal: pulse.env*cos(pulse.phi) (Tesla) 
%pulse.esd: Energy Spectral Density of pulse
%pulse.esdf: Frequency axis for ESD (Hz)
%pulse FWHM: FWHM of ESD (linear scale). This is the pulse bandwidth
%pulse.Qcrit: minimum adiabaticity of pulse (only for Adiabatic pulses)
%
%E.g.
%pulse.name='square';
%pulse.B1max=5.6e-4;
%pulse.Tp=32e-9;
%pulse.f_max=0;
%pulse=Create_Pulse(pulse);

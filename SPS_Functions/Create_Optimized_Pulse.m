% This function generates an optimized pulse for a user-specified
% operation.
% Syntax:
%[ pulse ] = Create_Optimized_Pulse( pulse )
% pulse is a structure with the following fields:
%
%--------------INPUTS (optional and vary depending on application):
% pulse.name: string. pulse name. e.g  'square', 'HSn', 'Hermite180'
% pulse.method: string. Optimization method. 'Fixed_B1', 'Fixed_Tp', 'Fixed_opBW',
% or 'None'
% pulse.regime: only for FM pulses. 'Adiabatic', or 'RapidPassage'
% pulse.npts: length of time vector
% pulse.phi_0: starting phase (radians)
% pulse.f_center: center frequency (Hz)
% pulse.Tp: pulse length (seconds) 
% pulse.B1max: B1 strength (Tesla)
% pulse.f_max: frequency sweep amplitude (Hz)
% pulse.opBW: target operation bandwidth (Hz) 
% pulse.theta: target rotation angle e.g. pi/2, pi. Invalid for pulses designed to perform specific
% operations (E.g. EBURP)
% pulse.R: Adiabatic design parameter (Time-bandwidth product). For Fixed_opBW optimization method
% pulse.K: Adiabatic design parameter for Fixed_B1 optimization method
% pulse.node: Adiabatic design parameter for chirp pulse 
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
%pulse.name='HSn'
%pulse.method='Fixed_opBW';
%pulse.regime='Adiabatic';
%pulse.theta=pi;
%pulse.opBW=60e6;
%pulse=Create_Optimized_Pulse(pulse);

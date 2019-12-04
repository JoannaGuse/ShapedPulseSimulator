
%This GUI exports pulse data to AWG compatible .wav files according to
%specified formatting inputs.
%
%examples:
%Export_wavfileGUI_Standalone(); % will result in a loading data window.
%Export_wavfileGUI_Standalone(pulse); %will load pulse 
%pulse must have the fields: pulse.Bx (I), pulse.By (Q), pulse.signal (Y),
%pulse.t

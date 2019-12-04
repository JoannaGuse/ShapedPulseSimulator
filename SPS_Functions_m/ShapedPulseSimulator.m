function [  ] = ShapedPulseSimulator(  )

%find 
str=which('ShapedPulseSimulator');
[filepath,~,~]=fileparts(str);
addpath(genpath(filepath))

disp('===========================================================')
disp('ShapedPulseSimulator is a toolbox for designing and simulating pulses for EPR spectropscopy.')
disp('Author: Joanna A. Guse')
disp('Version: ShapedPulseSimulator 1.1')
disp('Release Date: 18/12/2017')
disp(['Folder: ',filepath])
disp('===========================================================')
help SPS_Functions
help INBUILT
end


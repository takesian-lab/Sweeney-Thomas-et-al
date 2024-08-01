function [stim_units, V1_label, V2_label] = assign_stim(stimProtocol)  

%Assign stim labels to stim with >1 type of stim condition
if stimProtocol == 2 %Receptive Field
    stim_units = {'kHz', 'dB'};
    V1_label = 'Frequency (kHz)';
    V2_label = 'Intensity (dB)';
elseif stimProtocol == 3 %FM sweep
    stim_units = {'ov/s', 'dB'};
    V1_label = 'Rate (ov/s)';
    V2_label = 'Intensity (dB)';
elseif stimProtocol == 5 %SAM
    stim_units = {'Hz', ''};
    V1_label = 'Rate (Hz)';
    V2_label = 'Modulation Depth';
elseif stimProtocol == 6 %SAM freq
    stim_units = {'kHz', ''};
    V1_label = 'Frequency (kHz)';
    V2_label = 'Modulation Depth';
elseif stimProtocol == 12 %Spontaneous
    stim_units = [];
    V1_label = [];
    V2_label =[];
else
    error('New stim protocol')
end
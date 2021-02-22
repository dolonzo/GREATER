% eeg_checkburst () - This function is used to check for burst TMS delivery
% within a concurrent TMS-EEG dataset. Burst delivery is determined by any
% pulse frequency greater than 20 Hz. The first pulse in each burst is
% identified with an additional 'First' event in order to allow epoching of
% the entire burst.
%
% INPUTS:
%   EEG - EEG dataset containing TMS pulse events
%   pulse_event - cell array containing the event type coding TMS pulses
%
% OUTPUTS:
%   EEG - EEG dataset with the first pulse in each burst identified with
%          an additional event at the same latency with 'First' as the
%          event type. If the dataset does not contain bursts, the EEG
%          dataset is returned unchanged.
%   isburst - logical 0 or 1 if the dataset contains burst (1) or not (0)
%
% USAGE:
%   [EEG, isburst] = eeg_checkburst(EEG, {'1'})
%
% Author: David Sorensen, 2021

function [EEG, isburst] = eeg_checkburst(EEG, pulse_event)

fprintf('Checking for burst patterning...\n');

pulse_indices = find(ismember({EEG.event(:).type}, pulse_event));
pulse_latencies = [EEG.event(pulse_indices).latency];

pulse_lat_diffs = nan(1, length(pulse_latencies));
for i = 2:length(pulse_latencies)
    pulse_lat_diffs(i) = pulse_latencies(i) - pulse_latencies(i-1);
end

burst_cutoff = 45e-3*EEG.srate; %burst_cutoff set to 45 ms to account for jitter in pulse frequency
bool_burst_latency = pulse_lat_diffs<burst_cutoff;

if ~any(bool_burst_latency)
    isburst = 0;
    return;
end

isburst = 1;

firsts = find(~bool_burst_latency);
burst_count = zeros(1, length(firsts)-1);
for i = 1:length(firsts)-1
    burst_count(i) = firsts(i+1)-firsts(i);
end

if any(burst_count~=burst_count(1))
    warning('Warning: burst counts not consistent. Check for errant pulses in the dataset.');
end

for ind = 1:length(firsts)
    EEG.event(end+1) = EEG.event(pulse_indices(firsts(ind)));
    EEG.event(end).type = 'First';
end


end


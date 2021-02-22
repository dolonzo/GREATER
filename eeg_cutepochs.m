% eeg_cutepochs() - This function is designed to be used with the
% eeg_spliceepochs function. eeg_cutepochs removes epochs from continuous
% data and returns both the epoched data and the continuous dataset with
% epochs removed. In the continuous dataset, removed epochs are replaced by
% boundary events with a precut field identifying EpochCuts. These events
% are utilized later by the eeg_spliceepochs function.
%
% INPUTS:
%   EEGIN - EEG dataset (as an eeglab EEG structure) from which to remove epochs
%   eventtype - char array or cell array containing event type(s) for
%               timelocking the epochs
%   latencies - [start end] the start and end latencies, relative to the
%               timelocking event(s) in eventtype, of the epochs in seconds
%
% OUTPUTS:
%   EEGOUT - EEG dataset (as an eeglab EEG structure) with the epochs
%            removed
%   EpochEEG - EEG dataset containing the epoched EEG data that was removed
%
% USAGE:
%   [ContEEG, EpochEEG] = eeg_cutepochs(EEG, {'1'}, [-50, 145]);
%
% See Also:
%   pop_epoch; pop_rmdat; eeg_spliceepochs
%
% Author: David Sorensen, 2020

function [EEGOUT, EpochEEG] = eeg_cutepochs(EEGIN, eventtype, latencies)

%Identify the events that exist before we remove anything
if ~isfield(EEGIN.event, 'precut')
    [EEGIN.event(:).precut] = deal('1');
else
    emptyprecuts = cellfun(@isempty, {EEGIN.event(:).precut});
    [EEGIN.event(emptyprecuts).precut] = deal('1');
    previousEpochCuts = find(contains({EEGIN.event(:).precut}, 'EpochCut'));
    if ~isempty(previousEpochCuts)        
        insertfunc = @(a) insertBefore(a, 'EpochCut', 'Old');
        newEpochCutLabels = cellfun(insertfunc, {EEGIN.event(previousEpochCuts).precut}, 'UniformOutput', false);
        [EEGIN.event(previousEpochCuts).precut] = newEpochCutLabels{:};
    end
end

fprintf('Cutting epochs out of continuous dataset\n');

[EpochEEG, acceptedEvents] = pop_epoch(EEGIN, eventtype, latencies);
urepochnum = num2cell(1:length(EpochEEG.epoch));
[EpochEEG.epoch(:).urepochnum] = urepochnum{:};

all_eventtype = find(ismember({EEGIN.event(:).type}, eventtype));
acceptedEvents = all_eventtype(acceptedEvents);

%EEGOUT = pop_rmdat(EEGIN, eventtype, latencies, 1); Cannot just use rmdat:
%inexact limits
lat1 = round([EEGIN.event(acceptedEvents).latency]) + EpochEEG.xmin*EpochEEG.srate;
lat2 = round([EEGIN.event(acceptedEvents).latency]) + EpochEEG.xmax*EpochEEG.srate;
pnt_array = [lat1', lat2'];

EEGOUT = pop_select(EEGIN, 'nopoint', pnt_array);
%Find the boundary events that have empty precut fields
newbounds = find(strcmp('boundary', {EEGOUT.event(:).type}) & cellfun(@isempty, {EEGOUT.event(:).precut}));

%Check for overlapping epochs and throw an error
if length(newbounds) ~= size(EpochEEG.data, 3)
    error('The number of epochs extracted and the number of EpochCut events do not match. This can occur as a result of epochs with overlapping time windows.')
end
%This could also be rectified, but more difficult

%Label them as EpochCut events
[EEGOUT.event(newbounds).precut] = deal('EpochCut');
end

% Things that can break this code:--REVISIT THIS BEFORE PUBLISHING
% Overlapping epoch requests--Compare number of EpochCut events to number
% of epochs
% Epochs that would go over other boundary events--not included when using
% pop_epoch. Compare number of epochs to number of events of eventtype.
% Epochs extending beyond the dataset--not included when using pop_epoch.
% Compare number of epochs to number of events of eventtype.
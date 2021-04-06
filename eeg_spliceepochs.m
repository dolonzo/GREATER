% eeg_spliceepochs() - This function is designed to be used with the
% eeg_cutepochs function. eeg_sliceepochs puts previously epoched data back
% into a continuous dataset containing EpochCut events, returning the newly
% spliced dataset
%
% INPUTS:
%   EEGIN - EEG dataset (as an eeglab EEG structure) containing EpochCut
%           events from eeg_cutepochs; else, if passed alone, ALLEEG
%           structure containing both datasets from eeg_cutepochs
%   EpochEEG - EEG dataset containing epoched EEG data produced using
%              eeg_cutepochs
%
% OUTPUTS:
%   EEGOUT - EEG dataset (as an eeglab EEG structure) with the epochs
%            removed
%
% USAGE:
%   EEG = eeg_spliceepochs(ContEEG, EpochEEG);
%
% See also:
%   eeg_cutepochs; pop_epoch;
%
% Author: David Sorensen, 2020

function EEGOUT = eeg_spliceepochs(EEGIN, EpochEEG)

if nargin < 1
    help eeg_spliceepochs
    return
end

if nargin < 2
    geometry = {[1, 0.5], [1, 0.5]};
    uilist = {...
        {'Style', 'text', 'string', 'Index of first dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'first'}...
        {'Style', 'text', 'string', 'Index of second dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'second'} };
    [~, ~, ~, outstruct, ~] = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Splice epochs back into continuous data--eeg_spliceepochs()');
    if ~isempty(outstruct)
        EEGOUT = eeg_spliceepochs(INEEG1(str2num(outstruct.first)), INEEG1(str2num(outstruct.second)));
    end
    return
end

%Find all the EpochCut events
epochcuts = find(strcmp('EpochCut', {EEGIN.event(:).precut}));
if isempty(epochcuts)
    error('No EpochCut events found in the dataset.\nYou can only splice epochs back in to a dataset that were previously removed using eeg_cutepochs');
end
%Put the data from EpochEEG in place of the matching EpochCut events using
%the epoch numbering system
epochs_to_insert = [EpochEEG.epoch(:).urepochnum];
EEGOUT = EEGIN;
%EEGOUT.data = [EEGOUT.data, reshape(EpochEEG.data, EpochEEG.nbchan, EpochEEG.pnts*EpochEEG.trials)];
EEGOUT.data = [];
data = zeros(EEGOUT.nbchan, EEGOUT.pnts+EpochEEG.pnts*EpochEEG.trials);
latency = round(EEGOUT.event(epochcuts(epochs_to_insert(1))).latency);
data(:, 1:latency-1) = EEGIN.data(:, 1:latency-1);
fprintf('Splicing epochs back into continuous dataset\n');
for ind = 1:length(epochs_to_insert)
    latency = round(EEGOUT.event(epochcuts(epochs_to_insert(ind))).latency); %round because latencies are not forced to be integers
    duration = EEGOUT.event(epochcuts(epochs_to_insert(ind))).duration;
    %get the event information from the EEGOUT dataset
    affectedevents = [EEGOUT.event(:).latency] > latency;
    %Copy the data from the EpochEEG back in
    %EEGOUT.data(:, latency:latency+duration-1) = EpochEEG.data(:, :, ind);
    data(:, latency:latency+duration-1) = EpochEEG.data(:, :, ind);
    %Copy only the intervening data until we get to the last epoch to
    %splice back in
    if ind <= length(epochs_to_insert)-1
        nextlat = round(EEGOUT.event(epochcuts(epochs_to_insert(ind+1))).latency);
        gapdur = nextlat-latency;
        %EEGOUT.data(:, latency+duration:latency+duration+gapdur) = EEGIN.data(:, latency-(ind-1)*duration:nextlat-(ind-1)*duration);
        data(:, latency+duration:latency+duration+gapdur) = EEGIN.data(:, latency-(ind-1)*duration:nextlat-(ind-1)*duration);
    else
        enddur = size(EEGIN.data(:, latency-(ind-1)*duration:end), 2);
        %EEGOUT.data(:, latency+duration:latency+duration+enddur-1) = EEGIN.data(:, latency-(ind-1)*duration:end);
        data(:, latency+duration:latency+duration+enddur-1) = EEGIN.data(:, latency-(ind-1)*duration:end);
    end
    %change latency of the rest of the events
    caffected = num2cell([EEGOUT.event(affectedevents).latency] + duration);
    [EEGOUT.event(affectedevents).latency] = caffected{:};
    %Add back in the events from the EpochEEG
    epochevents = EpochEEG.epoch(ind).event;
    tmpevent = EpochEEG.event(epochevents);
    %Change the latencies of the events to match the latencies of the
    %continuous dataset
    ceventlatencies = num2cell(latency + eeg_lat2point([EpochEEG.epoch(ind).eventlatency], 1, EpochEEG.srate, [EpochEEG.xmin*1000 EpochEEG.xmax*1000], .001));
    [tmpevent.latency] = ceventlatencies{:};
    %Remove the epoch field from the event structure to make it compatible
    %with continuous data event structure
    tmpevent = rmfield(tmpevent, 'epoch');
    %urevent structure can cause problems, and does not make much sense in
    %this context
    if isfield(EEGOUT.event, 'urevent') && ~isfield(tmpevent, 'urevent')
        EEGOUT.event = rmfield(EEGOUT.event, 'urevent');
    elseif ~isfield(EEGOUT.event, 'urevent') && isfield(tmpevent, 'urevent')
        tmpevent = rmfield(tmpevent);
    end
    %Append the events to the end of the EEGOUT event structure
    EEGOUT.event = [EEGOUT.event, tmpevent];
end
%Delete the EpochCut events once they have been replaced
EEGOUT.event(epochcuts(epochs_to_insert)) = [];
notspliced = find(strcmp('EpochCut', {EEGOUT.event(:).precut}));
if ~isempty(notspliced)
    [EEGOUT.event(notspliced).precut] = deal('NonSplicedCut');
end
oldepochcuts = find(contains({EEGOUT.event(:).precut}, 'OldEpochCut'));
if ~isempty(oldepochcuts)
    extractfunc = @(a) extractAfter(a, 'Old');
    paredepochcuts = cellfun(extractfunc, {EEGOUT.event(oldepochcuts).precut}, 'UniformOutput', false);
    [EEGOUT.event(oldepochcuts).precut] = paredepochcuts{:};
end
EEGOUT.data = data;
EEGOUT.pnts = length(EEGOUT.data);
EEGOUT = eeg_checkset(EEGOUT);
EEGOUT = eeg_checkset(EEGOUT, 'eventconsistency');

end
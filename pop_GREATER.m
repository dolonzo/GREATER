% pop_GREATER() - Runs the GREATER analysis to remove large TMS pulse
% artifacts in a first round of spatial filtering. In brief, this pipeline
% consists of the following steps:
%       1. Remove and interpolate the pulse artifact
%       2. Downsample the data
%       3. Remove interpolated data prior to spatial filtering
%       4. Separate Intertrain Interval data from TMS pulse epochs
%       5. Reject bad channels based on FASTER criteria applied to the
%       Intertrain Interval data
%       6. Reject bad portions of data from Intertrain Interval data based
%       on FASTER criteria
%       7. Compute the Generalized Singular Value Decomposition (GSVD)
%       jointly on the Intetrain Interval data and the TMS Pulse epochs
%       8. (Optionally) Plot the components
%       9. (Optionally) Reject artifactual GSVD components from the
%       dataset. This can be done based on amplitude in a post-pulse
%       window or based on the singular values
%       10. Recombine the TMS pulse epochs with the Intertrain Interval
%       data for further processing and analysis
%
% INPUTS:
%   EEG - a TMS-EEG dataset as an eeglab EEG structure
%   pulse_event - the event type of the TMS pulse marker
%   pulse_window - [start end] a 2x1 numerical array containing the start 
%                   and end latency (in ms) of the pulse window to be
%                   removed and interpolated relative to the pulse event
%   prepulse - [start end] prepulse baseline (in ms) used to replace the
%               pulse window to avoid biasing the GSVD with interpolated
%               data
%   interp_window - [before after] length of data (in ms) before and after
%                    the pulse_window. Note that before should be given as
%                    a positive integer
%   sampling_rate - sampling rate (in Hz) to resample the data to
%   epoch_window - [start end] start and end times (in ms) relative to the
%                   pulse_event from which to extract epochs
%   ref_chan - the index of the reference channel for the data. If the
%               reference channel is not included in the data, pass an
%               index greater than the number of channels. An additional
%               channel will be added with zero at each time point.
%   savecomps - 0|1 save the GSVD components in the CONTEEG and EpochEEG
%                datasets. Defaults to 0 to save memory (components can be
%                computed when needed)
%   plotcomps - 0|1 plot the components after calculating
%   rej - 0|1 whether to reject and remove artifactual components
%   rej_criteria - 'amp'|'gsv' the type of criterion used to reject
%                    artifactual components. Either amplitude in the post
%                    pulse period or the singular values can be used
%   rej_thresh - numeric threshold (either amplitude or singular value)
%                 which leads to rejection of components
%   muscle_win - [duration] duration of time after the removed data over
%                 which to check the amplitude when using 'amp' rej_criteria
%
% OUTPUTS:
%   OUTEEG - the recombined EEG dataset with bad channels identified, bad
%             portions of Intertrain Interval data removed, and
%             (optionally) artifactual GSVD components removed. Because
%             GSVD is specific to the two datasets used, the GSVD-related
%             fields are removed from the dataset
%   CONTEEG - the Intertrain Interval data after the same processing
%              containing the GSVD weights and svs
%   EpochEEG - the Pulse epochs after the same processing containing the
%               GSVD weights and svs
%
% USAGE:
%   EEG = pop_GREATER(EEG); interactive popup window
%   [EEG, Intertrain, Pulse] = pop_GREATER(EEG, {'1'}, [-1, 10], [-25, -2],...
%                                [5, 5], 1024, [-50, 145], 64, 0, 1, 0);
%                                calculate the GSVD but do not reject
%                                components
%
%   EEG = pop_GREATER(EEG, {'1'}, [-1, 10], [-25, -2], [5, 5], 1024,...
%           [-50, 145], 64, 0, 1, 1, 'amp', 30, 20); reject components
%           based on post-pulse amplitude in the 20 ms after the tms pulse.
%           Components with amplitudes greater than 30 a.u. will be
%           rejected.
%
%   EEG = pop_GREATER(EEG, {'1'}, [-1, 10], [-25, -2], [5, 5], 1024,...
%           [-50, 145], 64, 0, 1, 1, 'gsv', 0.5); reject components based
%           on generalized singular values. Components with singular values
%           less than 0.5 (Intertrain to Pulse Epochs) will be removed.
%
% See Also:
%   pop_eeggsvd; pop_rejgsvdcomps_amp; pop_rejgsvdcomps_gsv;
%   pop_subgsvdcomps; eeg_cutepochs; eeg_spliceepochs
%
% DEPENDENCIES:
%   FASTER: channel_properties, min_z, epoch_properties
%   TESA: pop_tesa_removedata (tesa_removedata), pop_tesa_interpdata
%   (tesa_interpdata) We can package these together under a GNU General
%   Public License
%
% Authored by: David Sorensen, 2021

function [OUTEEG, CONTEEG, EpochEEG] = pop_GREATER(EEG, pulse_event, pulse_window, prepulse,  interp_window, sampling_rate, epoch_window, ref_chan, savecomps, plotcomps, rej, rej_criteria, rej_thresh, muscle_win)

if nargin < 1
    help pop_GREATER
    return
end

if nargin < 2
    geometry = { [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [1, 1, 2], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]};
    uilist = {...
        {'Style', 'text', 'string', 'Event(s) marking pulses'}...
        {'Style', 'edit', 'string', '' 'tag' 'event'} {}...
        {'Style', 'text', 'string', 'Window for removal [start end] (ms)'}...
        {'Style', 'edit', 'string', '' 'tag' 'pulsewin'} {}...
        {'Style', 'text', 'string', 'Prepulse baseline [start end] (ms)'}...
        {'Style', 'edit', 'string', '' 'tag' 'prepulse'} {}...
        {'Style', 'text', 'string', 'Window for interpolation [before after] (ms)'}...
        {'Style', 'edit', 'string', '' 'tag' 'interpwin'} {}...
        {'Style', 'text', 'string', 'New sampling rate'}...
        {'Style', 'edit', 'string', '' 'tag' 'srate'} {}...
        {'Style', 'text', 'string', 'Epoch window relative to pulse [start end] (ms)'}...
        {'Style', 'edit', 'string', '' 'tag' 'epochwin'} {}...
        {'Style', 'text', 'string', 'Reference channel index'}...
        {'Style', 'edit', 'string', '' 'tag' 'refchan'} {}...
        {'Style', 'checkbox', 'string', 'Save component activations' 'tag' 'save'} ...
        {'Style', 'checkbox', 'string', 'Plot components' 'tag' 'plot'} ...
        {'Style', 'checkbox', 'string', 'Reject components' 'tag' 'rej'}...
        {} {} {'Style', 'text', 'string', 'Rejection criterion'}...
        {'Style', 'popupmenu', 'string', {'Amplitude', 'Singular values'} 'tag' 'rejcrit'} ...
        {} {} {'Style', 'text', 'string', 'Threshold for rejection'}...
        {'Style', 'edit', 'string', '' 'tag' 'rejthresh'}...
        {} {} {'Style', 'text', 'string', 'Window for amplitude rejection'}...
        {'Style', 'edit', 'string', '' 'tag' 'rejwin'} };
    [~, ~, ~, outstruct, ~] = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'GREATER Pipeline--pop_GREATER()');
    if ~isempty(outstruct)
        rej_criterion = ['amp', 'gsv'];
        [OUTEEG, CONTEEG, EpochEEG] = pop_GREATER(EEG, cell(outstruct.event), str2num(outstruct.pulsewin), ...
            str2num(outstruct.prepulse), str2num(outstruct.interpwin), str2num(outstruct.srate), ...
            str2num(outstruct.epochwin), str2num(outstruct.refchan), outstruct.save, outstruct.plot, oustruct.rej, ...
            rej_criterion(outstruct.rejcrit), str2num(outstruct.rejthresh), str2num(outstruct.rejwin));
    end
    return
end
        
%Remove and interpolate data prior to downsampling
EEG = pop_tesa_removedata(EEG, pulse_window, prepulse, pulse_event);
EEG = pop_tesa_interpdata(EEG, 'cubic', interp_window);

%Downsample
EEG = pop_resample(EEG, sampling_rate);

%Check for burst epochs
[EEG, isburst] = eeg_checkburst(EEG, pulse_event);

%Add zero data reference channel if not included in the dataset
if ref_chan > EEG.nbchan
    fprintf("Adding reference channel 'REF' to dataset\n");
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.nbchan = EEG.nbchan+1;
    EEG = pop_chanedit(EEG, 'append', {EEG.nbchan 'REF' 0 0 0 0 0 0 0 0 [] [] [] 1});
    EEG = eeg_checkset(EEG);
    ref_chan = EEG.nbchan;
end

%Remove interpolated data
EEG = pop_tesa_removedata(EEG, pulse_window, prepulse, pulse_event);

%Separate out the pulse epoch data from intertrain interval data
if isburst
    burst_event = {'First'};
    [CONTEEG, EpochEEG] = eeg_cutepochs(EEG, burst_event, epoch_window./1000);
else
    [CONTEEG, EpochEEG] = eeg_cutepochs(EEG, pulse_event, epoch_window./1000);
end

%Select channels as bad based on intertrain interval data (using FASTER
%criteria)
list_properties = channel_properties(CONTEEG, 1:EEG.nbchan, ref_chan);
exceeded_threshold = min_z(list_properties);
gsvd_chans = setdiff(find(~exceeded_threshold), ref_chan);

%Remove bad segments from intertrain interval data (using a regular
%epoch window and FASTER epoch rejection criteria)
CONTEEG = eeg_regepochs(CONTEEG, 'extractepochs', 'off');
[tmpEEG, tmpITIEpochs] = eeg_cutepochs(CONTEEG, {'X'}, [-0.05, 0.9]);
list_properties = epoch_properties(tmpITIEpochs, gsvd_chans');
exceeded_threshold = min_z(list_properties);
tmpITIEpochs = pop_rejepoch(tmpITIEpochs, find(exceeded_threshold), 0);
CONTEEG = eeg_spliceepochs(tmpEEG, tmpITIEpochs);
clear tmpEEG tmpITIEpochs

%Calculate GSVD
[CONTEEG, EpochEEG] = pop_eeggsvd(CONTEEG, EpochEEG, savecomps, gsvd_chans);

%Plot GSVD components
if plotcomps
    pop_plotgsvdcomps(CONTEEG, 'title', 'Intertrain Interval Data GSVD Components');
    pop_plotgsvdcomps(EpochEEG, 'title', 'TMS Pulse Data GSVD Components');
    uiwait;
end

%Reject and remove components
if rej
    %Identify components to reject
    switch rej_criteria
        case 'amp'
            comps_to_sub = pop_rejgsvdcomps_amp(EpochEEG, rej_thresh, muscle_win);
        case 'gsv'
            comps_to_sub = pop_rejgsvdcomps_gsv(CONTEEG, EpochEEG, rej_thresh);
    end
    %Remove the component from the dataset
    [CONTEEG, EpochEEG] = pop_subgsvdcomps(CONTEEG, EpochEEG, comps_to_sub);
end

%Splice the data back together for subsequent processing and analysis
OUTEEG = eeg_spliceepochs(CONTEEG, EpochEEG);
OUTEEG = rmfield(OUTEEG, {'gsvdcomp', 'gsvdwts', 'gsvdsv'}); %keep gsvd_chans in to be used as channel set for further processing

end

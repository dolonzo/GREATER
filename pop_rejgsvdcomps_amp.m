% pop_rejgsvdcomps_amp() - function to mark GSVD components for rejection
%                          based on the amplitude in the post-pulse window
%                          (for TMS-EEG artifact removal) relative to the
%                          pre-pulse baseline.
%
% INPUTS:
%   PulseEpoch - EEG dataset containing TMS pulse epochs which have had the
%                major pulse artifact removed using pop_tesa_removedata
%   threshold - Threshold value above which to mark the component for
%               rejection
%   art_win - Time (in ms) relative to the end of the removed data over
%             which to average the amplitude
%
% OUTPUTS:
%   rej_comps - array containing component indices marked for rejection
%   crit_value - amplitude values for each component
%
% USAGE:
%   [comps_to_rej, amplitude] = pop_rejgsvdcomps_amp(EEG, 15, 20); rejects
%                   components whose average value in the 20 ms post pulse
%                   exceeds 15 arbitrary units
%
% See Also:
%   pop_eeggsvd; pop_rejgsvdcomps_gsv; pop_subgsvdcomps
%
% Authored by David Sorensen, 2021

function [rej_comps, crit_value] = pop_rejgsvdcomps_amp(PulseEpoch, threshold, art_win)

if nargin < 1
    help pop_rejgsvdcomps_amp
end

if nargin < 2
    geometry = {[1, 0.5], [1, 0.5]};
    uilist = {...
        {'Style', 'text', 'string', 'Rejection threshold'}...
        {'Style', 'edit', 'string', '' 'tag' 'thresh'}...
        {'Style', 'text', 'string', 'Post-pulse window endpoint (ms)'}...
        {'Style', 'edit', 'string', '', 'tag', 'artwin'} };
    [~, ~, ~, outstruct, ~] = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Post-pulse amplitude component rejection--pop_rejgsvdcomps_amp()');
    if ~isempty(outstruct)
        [rej_comps, crit_value] = pop_rejgsvdcomps_amp(PulseEpoch, str2num(outstruct.thresh), str2num(outstruct.artwin));
    end
    return
end

if ~isempty(PulseEpoch.gsvdcomp)
    gsvd_comp = PulseEpoch.gsvdcomp;
else
    gsvd_comp = PulseEpoch.data(PulseEpoch.gsvdchans, :)'*inv(PulseEpoch.gsvdwts')*inv(PulseEpoch.gsvdsv);
end

fprintf('Calculating amplitude in post-pulse window...\n');

componentERP = mean(reshape(diag(rms(PulseEpoch.gsvdwts))*(gsvd_comp*PulseEpoch.gsvdsv)', length(PulseEpoch.gsvdchans), PulseEpoch.pnts, PulseEpoch.trials), 3);
startremoved = eeg_lat2point(PulseEpoch.tmscut(end).cutTimesTMS(1)*1e-3, 1, PulseEpoch.srate, [PulseEpoch.xmin, PulseEpoch.xmax], 1); %check if cut epochs maintain TMS removed window
endremoved = eeg_lat2point(PulseEpoch.tmscut(end).cutTimesTMS(2)*1e-3, 1, PulseEpoch.srate, [PulseEpoch.xmin, PulseEpoch.xmax], 1);
removedpulse = mean(componentERP(:, round(startremoved):round(endremoved)), 2);
art_win_pnt = round(art_win*1e-3*PulseEpoch.srate);
artwindow = mean(componentERP(:, round(endremoved)+1:round(endremoved)+art_win_pnt), 2);

crit_value = abs(artwindow-removedpulse);
rej_comps = find(crit_value>threshold);

end


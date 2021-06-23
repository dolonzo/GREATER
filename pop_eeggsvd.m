% pop_eeggsvd() - This function is used to compute the GSVD between two
% EEGLAB EEG datasets. The two input datasets should contain the same
% number of channels.
%
% INPUTS:
%   INEEG1 - EEG dataset (as an eeglab EEG structure); alternatively, if
%            passed alone, the ALLEEG structure. A pop-up window will allow
%            selection of datasets within the ALLEEG structure.
%   INEEG2 - EEG dataset (as an eeglab EEG structure)
%   (optional) savecomps - [0|1] whether to save the GSVD components.
%                          Defaults to not saving components to reduce data
%                          size. Components can be recomputed later from
%                          data, weights, and sv matrices.
%   (optional) gsvd_chans - channel indices over which to compute the GSVD;
%                           defaults to all channels
%
% OUTPUTS:
%   OUTEEG1 - EEG dataset (as an eeglab EEG structure) with fields
%               containing the GSVD results
%   OUTEEG2 - EEG dataset (as an eeglab EEG structure) with fields
%               containing the GSVD results
% 
% USAGE:
%   [EEG1, EEG2] = pop_eeggsvd(EEG1, EEG2, 0, 1:63);
%
% See Also:
%   eeg_cutepochs; pop_GREATER
%
% Author: David Sorensen, 2021

function [OUTEEG1, OUTEEG2] = pop_eeggsvd(INEEG1, INEEG2, savecomps, gsvd_chans)

if nargin < 1
    help pop_eeggsvd
    return
end

%GUI window (if we can get to it)
if nargin == 1 %assume ALLEEG is passed
    uilist = {...
        {'Style', 'text', 'string', 'Index of first dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'first'}...
        {'Style', 'text', 'string', 'Index of second dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'second'}...
        {'Style', 'text', 'string', 'Channel indices to include in the GSVD'}...
        {'Style', 'edit', 'string', '', 'tag', 'chans'}... 
        {'Style', 'text', 'string', 'Save component activations?'}...
        {'Style', 'checkbox', 'string', '' 'tag' 'save'}};
    geometry = {[1, 0.5], [1, 0.5], [1, 0.5], [1, 0.5]};
    [~, ~, ~, outstruct, ~] = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Calculate GSVD--pop_eeggsvd()');
    if ~isempty(outstruct)
        [OUTEEG1, OUTEEG2] = pop_eeggsvd(INEEG1(str2num(outstruct.first)), INEEG1(str2num(outstruct.second)), outstruct.save, str2num(outstruct.chans));
    end
    return
end

if INEEG1.nbchan ~= INEEG2.nbchan
    error('The datasets must contain the same number of channels')
end

if nargin < 3
    savecomps = 0;
end

if nargin < 4
    gsvd_chans = 1:INEEG1.nbchan
end

INEEG1 = pop_rmbase(INEEG1, []);
INEEG2 = pop_rmbase(INEEG2, []);

fprintf('Calculating the GSVD...\n')

[comp1, comp2, weights, sv1, sv2] = gsvd(INEEG1.data(gsvd_chans, :)', INEEG2.data(gsvd_chans, :)', 0);

OUTEEG1 = INEEG1;
OUTEEG2 = INEEG2;

OUTEEG1.gsvdchans = gsvd_chans;
OUTEEG2.gsvdchans = gsvd_chans;

if savecomps
    OUTEEG1.gsvdcomp = comp1;
    OUTEEG2.gsvdcomp = comp2;
else
    OUTEEG1.gsvdcomp = [];
    OUTEEG2.gsvdcomp = [];
end

OUTEEG1.gsvdwts = weights;
OUTEEG2.gsvdwts = weights;
OUTEEG1.gsvdsv = sv1;
OUTEEG2.gsvdsv = sv2;
        
end

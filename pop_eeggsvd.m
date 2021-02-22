% pop_eeggsvd() - This function is used to compute the GSVD between two
% EEGLAB EEG datasets. The two input datasets should contain the same
% number of channels.
%
% INPUTS:
%   INEEG1 - EEG dataset (as an eeglab EEG structure)
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
%GUI window (if we can get to it)

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
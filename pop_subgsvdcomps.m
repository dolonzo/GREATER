% pop_subgsvdcomps() - This function subtracts the indicated GSVD
% components from both datasets used in calculating the decomposition.
%
% INPUTS:
%   INEEG1 - EEG dataset (as an eeglab EEG structure) used in calculating
%            the GSVD; else, passed alone, the ALLEEG structure containing
%            both EEG datasets used in calculating the GSVD
%   INEEG2 - EEG dataset (as an eeglab EEG structure) used in calculating
%            the GSVD;
%   comps_to_sub - array of component indices to be removed from the data
%
% OUTPUTS:
%   OUTEEG1 - EEG dataset from INEEG1 with the components removed
%   OUTEEG2 - EEG dataset from INEEG2 with the components removed
%
% USAGE:
%   [EEG1, EEG2] = pop_subgsvdcomps(ALLEEG); interactive pop up window
%   [EEG1, EEG2] = pop_subgsvdcomps(EEG1, EEG2, [1, 2, 4]); removes
%       components 1, 2, and 4 from both datasets
%
% See Also:
%   pop_eeggsvd; pop_rejgsvdcomps_gsv; pop_rejgsvdcomps_amp
%
% Authored by: David Sorensen, 2021

function [OUTEEG1,OUTEEG2] = pop_subgsvdcomps(INEEG1,INEEG2, comps_to_sub)

if nargin < 1
    help pop_subgsvdcomps
    return
end

%GUI window
if nargin == 1
    geometry = {[1, 0.5], [1, 0.5], [1, 0.5]};
    uilist = {...
        {'Style', 'text', 'string', 'Index of first dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'first'}...
        {'Style', 'text', 'string', 'Index of second dataset'}...
        {'Style', 'edit', 'string', '' 'tag' 'second'}...
        {'Style', 'text', 'string', 'Components to remove'}...
        {'Style', 'edit', 'string', '' 'tag' 'comps'} };
    [~, ~, ~, outstruct, ~] = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Subtract GSVD components--pop_subgsvdcomps()');
    if ~isempty(outstruct)
        [OUTEEG1, OUTEEG2] = pop_subgsvdcomps(INEEG1(str2num(outstruct.first)), INEEG1(str2num(outstruct.second)), str2num(outstruct.comps));
    end
    return
end

if isempty(INEEG1.gsvdcomp)
    gsvd_comp1 = INEEG1.data(INEEG1.gsvdchans, :)'*inv(INEEG1.gsvdwts')*inv(INEEG1.gsvdsv);
else
    gsvd_comp1 = INEEG1.gsvdcomp;
end

fprintf('Removing GSVD components ');
fprintf('%d ', comps_to_sub);
fprintf('from datasets\n');

S1 = INEEG1.gsvdsv;
S1(comps_to_sub, :) = 0;
OUTEEG1 = INEEG1;
OUTEEG1.data(INEEG1.gsvdchans, :) = (gsvd_comp1*S1*INEEG1.gsvdwts')';

if isempty(INEEG2.gsvdcomp)
    gsvd_comp2 = INEEG2.data(INEEG2.gsvdchans, :)'*inv(INEEG2.gsvdwts')*inv(INEEG2.gsvdsv);
else
    gsvd_comp2 = INEEG2.gsvdcomp;
end

S2 = INEEG2.gsvdsv;
S2(comps_to_sub, :) = 0;
OUTEEG2 = INEEG2;
OUTEEG2.data(INEEG2.gsvdchans, :) = (gsvd_comp2*S2*INEEG2.gsvdwts')';

end


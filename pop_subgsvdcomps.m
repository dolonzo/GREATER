% pop_subgsvdcomps() - This function subtracts the indicated GSVD
% components from both datasets used in calculating the decomposition.
%
% INPUTS:
%   INEEG1 - EEG dataset (as an eeglab EEG structure) used in calculating
%            the GSVD
%   INEEG2 - EEG dataset (as an eeglab EEG structure) used in calculating
%            the GSVD
%   comps_to_sub - array of component indices to be removed from the data
%
% OUTPUTS:
%   OUTEEG1 - EEG dataset from INEEG1 with the components removed
%   OUTEEG2 - EEG dataset from INEEG2 with the components removed
%
% USAGE:
%   [EEG1, EEG2] = pop_subgsvdcomps(EEG1, EEG2, [1, 2, 4]); removes
%       components 1, 2, and 4 from both datasets
%
% See Also:
%   pop_eeggsvd; pop_rejgsvdcomps_gsv; pop_rejgsvdcomps_amp
%
% Authored by: David Sorensen, 2021

function [OUTEEG1,OUTEEG2] = pop_subgsvdcomps(INEEG1,INEEG2, comps_to_sub)

if isempty(INEEG1.gsvdcomp)
    gsvd_comp1 = INEEG1.data(INEEG1.gsvdchans, :)'*inv(INEEG1.gsvdwts')*inv(INEEG1.gsvdsv);
else
    gsvd_comp1 = INEEG1.gsvdcomp;
end

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


% pop_rejgsvdcomps_gsv() - function to mark GSVD components for rejection
%                          based on the generalized singular values.
%
% INPUTS:
%   INEEG1 - EEG dataset used in calculating the GSVD using pop_eeggsvd
%   INEEG2 - EEG dataset used in calculating the GSVD using pop_eeggsvd
%   threshold - Threshold value relative to which the components will be
%               rejected. If the threshold is less than 1, components with
%               singular values smaller than the threshold will be marked
%               for rejection. If the threshold is greater than 1,
%               components with singular values greater than the threshold
%               will be marked for rejection.
%
% OUTPUTS:
%   rej_comps - array containing component indices marked for rejection
%
% USAGE:
%   comps_to_rej = pop_rejgsvdcomps_gsv(EEG1, EEG2, 0.5); returns indices
%           for components with generalized singular values less than 0.5
%   comps_to_rej = pop_rejgsvdcomps_gsv(EEG1, EEG2, 2); returns indices
%           for components with generalized singular values greater than 2
%
% See Also:
%   pop_rejgsvdcomps_gsv; pop_eeggsvd; pop_subgsvdcomps
% 
% Authored by David Sorensen, 2021

function rej_comps = pop_rejgsvdcomps_gsv(INEEG1,INEEG2, threshold)

% determine threshold direction relative to 1 (equal weighting)
if threshold > 1 
    gsv = diag(INEEG1.gsvdsv)./diag(INEEG2.gsvdsv);
    rej_comps = find(gsv>threshold);
else
    gsv = diag(INEEG1.gsvdsv)./diag(INEEG2.gsvdsv);
    rej_comps = find(gsv<threshold);
end

end


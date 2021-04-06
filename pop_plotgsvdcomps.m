% pop_plotgsvdcomps() - interactive plotting function for plotting GSVD
% components
%
% INPUTS:
%   EEG1 - EEG dataset containing computed GSVD components
%
%   eegplot options can be passed in 'key' 'val' pairs
%
% USAGE:
%   pop_plotgsvdcomps(EEG)
%   pop_plotgsvdcomps(EEG, 'winlength', 10, 'title', 'GSVD Components'...)
%
% See Also:
%   eegplot; pop_eeggsvd
%
% Authored by: David Sorensen, 2021

function pop_plotgsvdcomps(EEG1, varargin)

if ~isfield(EEG1, 'gsvdwts')
    error('No GSVD data are saved in this dataset');
end

if ~isempty(EEG1.gsvdcomp)
    gsvd_data = diag(rms(EEG1.gsvdwts))*(EEG1.gsvdcomp*EEG1.gsvdsv)';
else
    gsvdcomps = EEG1.data(EEG1.gsvdchans, :)'*inv(EEG1.gsvdwts')*inv(EEG1.gsvdsv);
    gsvd_data = diag(rms(EEG1.gsvdwts))*(gsvdcomps*EEG1.gsvdsv)';
end

eegplot(gsvd_data, 'srate', EEG1.srate, 'events', EEG1.event, varargin{:}, 'title', 'GSVD components');

end


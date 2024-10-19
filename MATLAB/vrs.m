function [vrs, f_n_out] = vrs(spectrum, f_n)
% vrms(spectrum) plots the vibration response spectrum of the base
% excitation given in spectrum. The first column of spectrum is an array of
% breakpoint frequencies, the second column of spectrum is an array of PSD
% values.
%
% [vrs, f_n_out] = vrs(spectrum, f_n) outputs the grms levels for each
% natural frequency in f_n_out. if f_n is not specified, f_n_out is 500
% logarithmically spaced points from the minimum to maximum frequency in
% spectrum. if f_n is supplied, f_n_out = f_n.

if nargin == 1
    f_n = logspace(log10(min(spectrum(:,1))), log10(max(spectrum(:,1))), 500)';
end

f_q = f_n';

sdof_response = sdof_PSD_response(spectrum, f_n, 10, f_q);
grms_values = grms(sdof_response, f_q);

if nargout == 0
    loglog(f_n, grms_values)
    grid on
    xlabel("Natural Frequency (Hz)")
    ylabel("Accel (GRMS)")
    title("Vibration Response Spectrum, SDOF Systems, Q=10")
    ylim([min(grms_values), 2^nextpow2(max(grms_values))])
else
    vrs = grms_values;
    f_n_out = f_n;
end
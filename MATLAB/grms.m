function grms = grms(spectrum, freqs)
% grms = GRMS(spectrum) calculates the grms value of a acceleration power spectrum
% where the first column of spectrum is an array of breakpoint frequencies,
% and the second column of spectrum is an array of APSD values.
%
% grms = GRMS(spectrum, freqs) returns a column vector of grms values given by the
% APSD values in each row of spectrum. freqs can be an array whose length =
% width(spectrum) or a matrix of size = size(spectrum).if freqs is an
% array, the ith element of freqs corresponds to the ith column of
% spectrum. if freqs is a matrix, the ith element of freqs corresponds to
% the ith element of spectrum (i.e. each row of spectrum has its own
% frequency array).

switch nargin
    case 1
        area = spectrum_area(spectrum);
    case 2
        area = spectrum_area(spectrum, freqs);
    otherwise
        error("grms takes 1 or 2 input arguments")
end

grms = sqrt(area);
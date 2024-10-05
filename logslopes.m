function slopes = logslopes(spectrum, freqs)
% slopes = logslopes(spectrum) returns an array of slopes between each 
% coordinate in a loglog power spectrum. the first column of spectrum is an
% array of breakpoint frequencies, and the second column of spectrum is an
% array of breakpoint PSD values.
%
% slopes = logslopes(spectrum, freqs) returns an matrix of slopes between
% each coordinate in the log-log power spectrum given by each row of 
% spectrum and the corresponding frequencies in freqs. freqs can be an 
% array whose length = width(spectrum) or a matrix of 
% size = size(spectrum). If freqs is an array, the ith element of freqs 
% corresponds to the ith column of spectrum. if freqs is a matrix, the ith 
% element of freqs corresponds to the ith element of spectrum (i.e. each 
% row of spectrum has its own frequency array).

switch nargin
    case 1
        frequency_consecutive_ratios = spectrum(2:end, 1)./spectrum(1:end-1, 1);
        g_square_consecutive_ratios = spectrum(2:end, 2)./spectrum(1:end-1, 2);
    case 2
        if size(freqs, 2) == 1
            freqs = freqs';
        end
        g_square_consecutive_ratios = spectrum(:, 2:end)./spectrum(:, 1:end-1);
        frequency_consecutive_ratios = freqs(:, 2:end)./freqs(:, 1:end-1);
    otherwise
        error("logslopes requires 1 or 2 input arguments")
end

slopes = log10(g_square_consecutive_ratios)./log10(frequency_consecutive_ratios);



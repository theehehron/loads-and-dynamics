function y_vals = spectrum_points(spectrum, freq_query)
% y_vals = spectrum_points(spectrum, freq_query) calculates y coordinates 
% on a loglog power spectrum corresponding to the frequencies in 
% freq_query. The first column of spectrum is the breakpoint frequencies,
% and the second column in spectrum is an array of breakpoint PSD values.

arguments
    spectrum (:, 2)
    freq_query double
end

slopes = logslopes(spectrum);

segment_equations = @(f, index) spectrum(index, 2)/(spectrum(index, 1)^slopes(index))*f^slopes(index);

y_vals = NaN(size(freq_query));
for i = 1:numel(freq_query)
    for j = 1:height(spectrum)-1
        if round(freq_query(i), 10) <= spectrum(j+1, 1) && freq_query(i) >= spectrum(j, 1)
            y_vals(i) = segment_equations(freq_query(i), j);
        end
    end
end
function sdof_response = sdof_PSD_response(base_spectrum, f_n, Q, f_query)
% SDOF_PSD_RESPONSE(base_spectrum, f_n, Q, f_query) given a base excitation 
% spectrum base_spectrum, sdof_response computes the PSD response of SDOF 
% systems with natural frequencies given by f_n at the frequencies given in
% f_query. The output is a length(f_n)xlength(f_query) matrix, where the 
% ith row is the response of an sdof oscilator with the ith natural 
% frequency in f_n. Sounds complicated but its not that bad. f_query 
% may also be a matrix, in which case, each row of sdof_response
% corresponds to the frequencies in each row of f_query.

if size(f_n, 2) ~= 1
    f_n = f_n';
end

[m, n] = size(f_query);
if m > 1 && n == 1
    f_query = f_query';
end

rho = f_query./f_n;
zeta = 1/(2*Q);

base_spectrum_y_values = spectrum_points(base_spectrum, f_query);

PSD_transfer_function = (1+(2.*zeta.*rho).^2) ./ ((1-rho.^2).^2+(2*zeta.*rho).^2);

sdof_response = PSD_transfer_function.*base_spectrum_y_values;
function area = spectrum_area(varargin)
% area = SPECTRUM_AREA(spectrum) calculates the area under an APSD curve given by 
% spectrum, where the first column of spectrum is an array of breakpoint 
% frequencies, and the second column of spectrum is an array of APSD 
% values.
%
% area = SPECTRUM_AREA(spectrum, freqs) returns a column vector of areas 
% under the APSD curves given in each row of spectrum and the corresponding
% frequencies in freqs. freqs can be an array whose length = 
% width(spectrum) or a matrix of size = size(spectrum). If freqs is an
% array, the ith element of freqs corresponds to the ith column of
% spectrum. if freqs is a matrix, the ith element of freqs corresponds to
% the ith element of spectrum (i.e. each row of spectrum has its own
% frequency array).

switch nargin
    case 1
        spectrum = varargin{1};
        n = logslopes(spectrum)';
        y = spectrum(:, 2)';
        f = spectrum(:, 1)';
    case 2
        y = varargin{1};
        f = varargin{2};
        if size(f, 2) == 1
            f = f';
        end
        n = logslopes(y, f);
    otherwise
        error("spectrum_area takes 1 or 2 input arguments")
end

segmentareas = zeros(size(n));
for i = 1:width(n)
    if n(i) == -1
        segmentareas(:,i) = y(:,i).*f(:,i).*log(f(:,i+1)./f(:,i));
    else
        segmentareas(:,i) = y(:,i)./(f(:,i).^n(:,i)).*(1./(n(:,i)+1)).*(f(:,i+1).^(n(:,i)+1) - f(:,i).^(n(:,i)+1));
    end
end

area = sum(segmentareas, 2);

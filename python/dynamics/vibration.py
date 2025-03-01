import numpy as np
from scipy import signal, fft

def loglogslopes(*args):
    """computes the slope between points on a loglog scale"""
    if len(args) == 1:
        spectrum = args[0]
        freq_consecutive_ratios = spectrum[1:,0]/spectrum[0:-1, 0]
        g_consecutive_ratios = spectrum[1:,1]/spectrum[0:-1, 1]
    elif len(args) == 2:
        y = np.atleast_2d(args[0])
        x = np.atleast_2d(args[1])
        g_consecutive_ratios = y[:, 1:]/y[:, 0:-1]
        freq_consecutive_ratios = x[:, 1:]/x[:, 0:-1]
    else:
        raise ValueError("loglogslopes takes one or two input arguments.")
    slopes = np.log10(g_consecutive_ratios)/np.log10(freq_consecutive_ratios)
    return slopes



def spectrum_area(*args):
    """calculates the area underneath a loglog plot"""
    if len(args) == 1:
        spectrum = args[0]
        n = np.atleast_2d(loglogslopes(spectrum))
        y = np.atleast_2d(spectrum[:,1])
        f = np.atleast_2d(spectrum[:,0])
    elif len(args) == 2:
        y = np.atleast_2d(args[0])
        f = np.atleast_2d(args[1])
        n = np.atleast_2d(loglogslopes(y, f))

    segment_areas = np.zeros(y.shape)
    
    if f.shape[0] == 1:
        f = f*np.ones((y.shape[0], 1))
    
    n_neg_one_mask = np.hstack((n == -1, np.zeros((n.shape[0], 1), dtype=bool)))
    n_non_neg_one_mask = np.hstack((n != -1, np.zeros((n.shape[0], 1), dtype=bool)))
    
    n_neg_one_mask_i_plus_1 = np.hstack((np.zeros((n.shape[0], 1), dtype=bool), n_neg_one_mask[:, :-1]))
    n_non_neg_one_mask_i_plus_1 = np.hstack((np.zeros((n.shape[0], 1), dtype=bool), n_non_neg_one_mask[:, :-1]))

    segment_areas[n_neg_one_mask] = y[n_neg_one_mask]*f[n_neg_one_mask]*np.log(f[n_neg_one_mask_i_plus_1]/f[n_neg_one_mask])
    segment_areas[n_non_neg_one_mask] = y[n_non_neg_one_mask]/(np.power(f[n_non_neg_one_mask], n[n!=-1]))*(1/(n[n!=-1]+1))* \
                            ((np.power(f[n_non_neg_one_mask_i_plus_1], (n[n!=-1]+1)))-(np.power(f[n_non_neg_one_mask], (n[n!=-1]+1))))
    
    area = np.sum(segment_areas, axis=1)
    return area



def grms(*args):
    """calculates the grms level of an acceleration power spectral density"""
    if len(args) == 1:
        area = spectrum_area(args[0])
    elif len(args) == 2:
        area = spectrum_area(args[0], args[1])
    else:
        raise ValueError("grms takes one or two input arguments.")
    
    grms = np.sqrt(area)
    return grms
    


def findnearest_above(x_q, x):
    """finds the element in x that is nearest but greater than x_q"""
    indices = np.searchsorted(x, x_q, side='right')
    indices = np.clip(indices, 0, len(x) - 1)
    
    return indices



def findnearest_below(x_q, x):
    """finds the element in x that is nearest and less than x_q"""
    indices = np.searchsorted(x, x_q, side='left')-1
    indices = np.clip(indices, 0, len(x) - 1)
    
    return indices



def spectrum_points(spectrum, freq_query):
    """given a power spectral density defined by the breakpoints in spectrum, 
    interpolates the spectrum y values at the x values in freq_query"""
    ## This works but is failing at endpoints. Need to fix.
    spectrum_freqs = spectrum[:,0]
    spectrum_accel = spectrum[:,1]
    freq_query = np.clip(freq_query, min(spectrum_freqs), max(spectrum_freqs))
    nearest_freq_above = spectrum_freqs[findnearest_above(freq_query, spectrum_freqs)]
    nearest_freq_below = spectrum_freqs[findnearest_below(freq_query, spectrum_freqs)]
    nearest_accel_above = spectrum_accel[findnearest_above(freq_query, spectrum_freqs)]
    nearest_accel_below = spectrum_accel[findnearest_below(freq_query, spectrum_freqs)]
    
    slopes = np.log10(nearest_accel_above/nearest_accel_below)/np.log10(nearest_freq_above/nearest_freq_below)
    spectrum_points = nearest_accel_below/np.power(nearest_freq_below, slopes)*np.power(freq_query, slopes)
    return spectrum_points
    


def sdof_psd_response(base_spectrum, f_n, Q, f_query):
    """Calculates the power spectrum of the response of an array of SDOF
    systems with natural frequencies f_n to the base excitation of base_spectrum"""
    f_n = np.atleast_2d(f_n)
    f_query = np.atleast_2d(f_query)
    
    if f_n.shape[1] != 1:
        f_n = f_n.T
    
    # if f_query.shape[0] > 1 & f_query.shape[1] == 1:
    #     f_query = f_query.T
    
    rho = f_query/f_n
    zeta = 1/(2*Q)
    
    base_spectrum_y_vals = spectrum_points(base_spectrum, f_query)
    psd_transfer_function = (1+np.power(2*zeta*rho, 2)) / ( np.power((1-np.power(rho, 2)), 2) + np.power(2*zeta*rho, 2) )
    
    sdof_response = psd_transfer_function*base_spectrum_y_vals
    return sdof_response
    

def vrs(*args):
    """calculates the vibration response spectrum to the input base excitation spectrum"""
    spectrum = args[0]
    if len(args) == 1:
        f_n = np.logspace(np.log10(min(spectrum[:,0])), np.log10(max(spectrum[:,0])), 500)
        f_n[-1] = spectrum[-1, 0] # set last element to equal largest element in spectrum
    elif len(args) == 2:
        f_n = args[1]
    else:
        raise ValueError("vrs takes one or two input arguments.")
        
    f_q = f_n.T 
    
    sdof_response = sdof_psd_response(spectrum, f_n, 10, f_q)
    vrs = grms(sdof_response, f_q);
    return vrs, f_q
    
    
    
def vrs_shock_equivalent(*args):
    """calculates the n-sigma vibration response spectrum to the input base excitation spectrum
    This is an 'equivalent shock response spectrum' of a random vibe signal described by the input spectrum"""
    spectrum = args[0]
    if len(args) == 1:
        f_n = np.logspace(np.log10(min(spectrum[:,0])), np.log10(max(spectrum[:,0])), 500)
        f_n[-1] = spectrum[-1, 0] # set last element to equal largest element in spectrum
        T = 60
    elif len(args) == 2:
        T = args[1]
        f_n = np.logspace(np.log10(min(spectrum[:,0])), np.log10(max(spectrum[:,0])), 500)
        f_n[-1] = spectrum[-1, 0] # set last element to equal largest element in spectrum
    elif len(args) == 3:
        T = args[1]
        f_n = args[2]
    else:
        raise ValueError("vrs_srs_equivalent takes one, two, or three input arguments.")
        
    f_q = f_n.T 
    
    nsigma = np.sqrt(2*np.log(f_n*T))
    sdof_response = sdof_psd_response(spectrum, f_n, 10, f_q)
    vrs = grms(sdof_response, f_q);
    vrs_nsigma = vrs*nsigma
    return vrs_nsigma, f_q



def vrs_miles(spectrum, Q, duration):
    """calculates estimated VRS using miles equation approximation, reference NASA SMC-S-016"""
    f_n = np.logspace(np.log10(min(spectrum[:,0])), np.log10(max(spectrum[:,0])), 500)
    f_n[-1] = spectrum[-1, 0] # set last element to equal largest element in spectrum
    
    G = spectrum_points(spectrum, f_n)
    n = np.sqrt(2*np.log(f_n*duration))
    vrs = n*miles(f_n, G, Q)
    return vrs, f_n
    

    
def miles(f_n, G, Q):
    grms = np.sqrt(np.pi/2*G*f_n*Q)
    return grms



def synthesize_vibration(spectrum, duration):
    """Synthesizes a random vibration time history described by the power spectral density in spectrum"""
    # 1. Generate low passed white noise time history.
    fs = 4000 # (Hz)
    T = duration # (sec)
    N = int(fs*T)

    mean = 0
    std_dev = 19
    white_noise = np.random.normal(mean, std_dev, N)
    
    # 2. Lowpass the white noise
    lowpass_filt_2000Hz = signal.butter(10, fs/2-1, 'lp', fs=fs, output='sos')
    white_noise_filtered = signal.sosfilt(lowpass_filt_2000Hz, white_noise)
    time = np.arange(N)/fs


    # 3. Compute the fft of the white noise time history.
    white_fft = fft.fft(white_noise_filtered)
    fft_freqs = fft.fftfreq(N, 1/fs)


    # 4. Calculate the target psd values that the white noise needs to be scaled to.
    PSD_target = spectrum_points(spectrum, np.abs(fft_freqs))


    # 5. scale the white noise signal in the frequency domain to match the target PSD
    scale_fft = np.sqrt(PSD_target*fs/(2*np.power(std_dev, 2)))
    synthesized_fft = scale_fft*white_fft


    # 6. Compute the inverse fft to get the time history
    synthesized_time_history = np.real(fft.ifft(synthesized_fft))
    
    return synthesized_time_history, time
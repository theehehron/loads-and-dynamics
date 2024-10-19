import numpy as np

def loglogslopes(*args):
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
    if len(args) == 1:
        area = spectrum_area(args[0])
    elif len(args) == 2:
        area = spectrum_area(args[0], args[1])
    else:
        raise ValueError("grms takes one or two input arguments.")
    
    grms = np.sqrt(area)
    return grms
    


def findnearest_above(x_q, x):

    indices = np.searchsorted(x, x_q, side='right')
    indices = np.clip(indices, 0, len(x) - 1)
    
    return indices



def findnearest_below(x_q, x):

    indices = np.searchsorted(x, x_q, side='left')-1
    indices = np.clip(indices, 0, len(x) - 1)
    
    return indices



def spectrum_points(spectrum, freq_query):
    ## This works but is failing at endpoints. Need to fix.
    spectrum_freqs = spectrum[:,0]
    spectrum_accel = spectrum[:,1]
    nearest_freq_above = spectrum_freqs[findnearest_above(freq_query, spectrum_freqs)]
    nearest_freq_below = spectrum_freqs[findnearest_below(freq_query, spectrum_freqs)]
    nearest_accel_above = spectrum_accel[findnearest_above(freq_query, spectrum_freqs)]
    nearest_accel_below = spectrum_accel[findnearest_below(freq_query, spectrum_freqs)]
    
    slopes = np.log10(nearest_accel_above/nearest_accel_below)/np.log10(nearest_freq_above/nearest_freq_below)
    spectrum_points = nearest_accel_below/np.power(nearest_freq_below, slopes)*np.power(freq_query, slopes)
    return spectrum_points
    


def sdof_psd_response(base_spectrum, f_n, Q, f_query):
    f_n = np.atleast_2d(f_n)
    f_query = np.atleast_2d(f_query)
    
    if f_n.shape[1] != 1:
        f_n = f_n.T
    
    if f_query.shape[0] > 1 & f_query.shape[1] == 1:
        f_query = f_query.T
    
    rho = f_query/f_n
    zeta = 1/(2*Q)
    
    base_spectrum_y_vals = spectrum_points(base_spectrum, f_query)
    psd_transfer_function = (1+np.power(2*zeta*rho, 2)) / ( np.power((1-np.power(rho, 2)), 2) + np.power(2*zeta*rho, 2) )
    
    sdof_response = psd_transfer_function*base_spectrum_y_vals
    return sdof_response



def vrs(*args):
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
    
    
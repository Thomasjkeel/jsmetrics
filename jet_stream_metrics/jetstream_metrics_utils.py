# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to identify or classify jet-stream in the literature
"""

### imports
import numpy as np
import xarray as xr
### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_resultant_wind(u_data, v_data):
    """
        TODO: work out where and how this will be used
    """
    return np.sqrt( u_data**2 + v_data**2 )


def get_latitude_and_speed_where_max_ws(data_row, latitude_col='lat'):
    """
        Write function description
    """
    if not data_row.isnull().all():
        max_speed_loc = data_row.argmax().data
        max_speed = data_row[max_speed_loc]
        lat_at_max = float(max_speed['lat'].values)
        speed_at_max = float(max_speed.data)
        return lat_at_max, speed_at_max 
    else:
        return None, None


def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    
    A low-pass filter removes short-term random fluctations in a time series

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    TAKEN FROM: https://scitools.org.uk/iris/docs/v1.2/examples/graphics/SOI_filtering.html

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[0+(window%2):-1] # edited from w[1:-1]


    from scipy import fftpack


def fourier_filter(data, timestep=1):
    """
        Carries out a Fourier transform for high frequency filtering
        TAKEN FROM: https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html
        NOTE: NOT CURRENTLY WORKING PROPERLY
        
        Parameters
        ----------
        data : np.array (1-d) 
            time series data
        timestep : float or int
            number used in the Discrete Fourier Transform sample frequencies (fftfreq)
    """
    fourier_transform = fftpack.fft(data)
    
    # The corresponding frequencies TODO: what does this do?
    sample_freq = fftpack.fftfreq(data.size, d=timestep)
    
    # And the power (sig_fft is of complex dtype)
    power = np.abs(data)**2
    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]
    
    high_freq_fft = fourier_transform.copy()
    high_freq_fft[np.abs(sample_freq) > peak_freq] = 0
    filtered_sig = fftpack.ifft(high_freq_fft)
    return filtered_sig
    
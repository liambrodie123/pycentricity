"""

A file for maximising overlaps for the pySEOBNRe package.

"""

import numpy as np
import bilby.gw.utils as utils
import bilby as bb
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy import signal
from scipy.optimize import minimize
import random


def plot_2d_overlap(overlaps, time_grid_mesh, phase_grid_mesh):
    """
    Plot the 2D grid of time vs phase shifts.
    :param overlaps: 2D grid
        grid of calculated overlaps
    :param time_grid_mesh: 1D grid
        grid of proposed time shift
    :param phase_grid_mesh: 1D grid
        grid of proposed phase shifts
    """
    plt.contourf(phase_grid_mesh, time_grid_mesh, overlaps)
    plt.ylabel("Time shift")
    plt.xlabel("Phase shift")
    plt.colorbar()
    plt.title("Overlap")
    plt.show()


def fourier_transform(waveform_time_domain, sampling_frequency):
    """
    Function for performing a Fourier transform on a dictionary of waveform polarisations.
    :param waveform_time_domain: dict or list of dicts
        single or multiple time-domain waveform polarisations
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :return:
        waveform_frequency_domain: dict
            frequency-domain waveform polarisations
    """
    if type(waveform_time_domain) is list:
        waveform_frequency_domain = []
        for single_waveform in waveform_time_domain:
            waveform_frequency_domain.append(
                {
                    key: bb.core.utils.nfft(
                        single_waveform[key], sampling_frequency=sampling_frequency
                    )[0]
                    for key in ["plus", "cross"]
                }
            )
    else:
        waveform_frequency_domain = {
            key: bb.core.utils.nfft(
                waveform_time_domain[key], sampling_frequency=sampling_frequency
            )[0]
            for key in ["plus", "cross"]
        }
    return waveform_frequency_domain


def wrap_by_n_indices(shift, waveform):
    """
    Wrap a time-domain waveform by a certain number of indices.
    :param shift: int
        the number of indices to shift the waveform by
    :param waveform: dict
        time-domain waveform polarisations
    :return:
        new_waveform: dict
            shifted time-domain waveform polarisations
    """
    new_waveform = deepcopy(waveform)
    new_waveform["plus"] = np.roll(new_waveform["plus"], shift=shift)
    new_waveform["cross"] = np.roll(new_waveform["cross"], shift=shift)
    return new_waveform


def wrap_at_maximum(waveform, max_index_other_model):
    """
    Wrap a time-domain waveform, shifting its maximum to the same point as the maximum of
    some comparison waveform.
    :param waveform: dict
        time-domain waveform polarisations
    :param max_index_other_model: int
        index of the maximum point of the comparison waveform
    :return:
        waveform: dict
            time-domain waveform polarisations
        shift: int
            the number of indices the waveform has been shifted by
    """
    max_index = np.argmax(np.abs(np.abs(waveform["plus"]) + np.abs(waveform["cross"])))
    shift = max_index_other_model - max_index
    waveform = wrap_by_n_indices(shift=shift, waveform=waveform)
    return waveform, shift


def overlap_function(
    a, b, frequency, psd, minimum_frequency=20, maximum_frequency=1024
):
    """
    Calculate the overlap between two waveforms.
    :param a: dict
        frequency-domain waveform polarisations of one waveform
    :param b: dict
        frequency-domain waveform polarisations of the other waveform
    :param frequency: array
        frequency array
    :param psd: array
        power spectral density array
    :param minimum_frequency: int
        minimum frequency at which to perform the calculation
    :param maximum_frequency: int
        maximum frequency at which to perform the calculation
    :return:
        overlap: float
            value of the overlap between the two waveforms
    """
    psd_interp = psd.power_spectral_density_interpolated(frequency)
    duration = 1.0 / (frequency[1] - frequency[0])
    minimum_frequency_index = np.where(frequency >= minimum_frequency)[0][0]
    maximum_frequency_index = np.where(frequency >= maximum_frequency)[0][0]
    # Defining temporary arrays to use
    _psd_interp = psd_interp[minimum_frequency_index:maximum_frequency_index]
    _a = {
        key: a[key][minimum_frequency_index:maximum_frequency_index] for key in a.keys()
    }
    _b = {
        key: b[key][minimum_frequency_index:maximum_frequency_index] for key in b.keys()
    }
    # Doing the calculation
    inner_a = utils.noise_weighted_inner_product(
        _a["plus"], _a["plus"], _psd_interp, duration
    )
    inner_a += utils.noise_weighted_inner_product(
        _a["cross"], _a["cross"], _psd_interp, duration
    )
    inner_b = utils.noise_weighted_inner_product(
        _b["plus"], _b["plus"], _psd_interp, duration
    )
    inner_b += utils.noise_weighted_inner_product(
        _b["cross"], _b["cross"], _psd_interp, duration
    )
    inner_ab = utils.noise_weighted_inner_product(
        _a["plus"], _b["plus"], _psd_interp, duration
    )
    inner_ab += utils.noise_weighted_inner_product(
        _a["cross"], _b["cross"], _psd_interp, duration
    )
    overlap = inner_ab / np.sqrt(inner_a * inner_b)
    return overlap.real


def zero_pad_frequency_domain_signal(waveform_frequency_domain, interferometers):
    """
    Mitigate the effects of spectral leakage
    :param waveform_frequency_domain: dict
        frequency-domain waveform polarisations
    :param interferometers: InterferometerList
        list of interferometers involved in the detection
    :return:
        waveform_frequency_domain: dict:
            frequency-domain waveform polarisations
    """
    ifo = interferometers[0]
    indices_to_zero_pad_low = np.where(ifo.frequency_array < ifo.minimum_frequency)[0]
    indices_to_zero_pad_high = np.where(ifo.frequency_array > ifo.maximum_frequency)[0]
    for key in waveform_frequency_domain.keys():
        waveform_frequency_domain[key][indices_to_zero_pad_low] = 0
        waveform_frequency_domain[key][indices_to_zero_pad_high] = 0
    return waveform_frequency_domain


def process_signal(waveform, comparison_length):
    """
    Zero-pad a waveform to a certain length, or truncate it if it is too long.
    :param waveform: dict
        time-domain waveform polarisations
    :param comparison_length: int
        number of indices we wish our waveform to have
    :return:
        waveform: dict
            zero-padded time-domain waveform polarisations
    """
    waveform_length = len(waveform['plus'])
    if waveform_length != comparison_length:
        if waveform_length < comparison_length:
            n_indices_to_pad = comparison_length - waveform_length
            waveform = {
                key: np.concatenate(([0] * n_indices_to_pad, waveform["plus"])) for key in waveform.keys()
            }
        elif waveform_length > comparison_length:
            n_indices_to_remove = waveform_length - comparison_length
            waveform  = {
                key: waveform[key][n_indices_to_remove:] for key in waveform.keys()
            }
    return waveform


def apply_tukey_window(waveform, alpha):
    """
    Return a time-domain waveform that has had a Tukey window applied to it.
    :param waveform: dict
        time-domain waveform polarisations
    :param alpha: float, 0 <= alpha <= 1
        rise gradient of the Tukey window
    :return:
        waveform: dict
            time-domain waveform polarisations
    """
    # First add some zeroes onto the upper end of the time series, to avoid the merge being swallowed by the window
    waveform = {
        key: np.concatenate((waveform[key], [0] * int(len(waveform['plus']) * alpha))) for key in waveform.keys()
    }
    waveform_length = len(waveform['plus'])
    waveform = dict(
        plus=np.multiply(
            signal.tukey(M=waveform_length, alpha=alpha), waveform["plus"]
        ),
        cross=np.multiply(
            signal.tukey(M=waveform_length, alpha=alpha), waveform["cross"]
        ),
    )
    return waveform


def apply_shifts(waveform_time_domain, index_shift, phase_shift, sampling_frequency):
    """
    Apply a given index shift and phase shift to a signal.
    :param waveform_time_domain: dict
        time-domain waveform polarisations
    :param index_shift: int
        number of indices to shift the waveform by
    :param phase_shift: float
        phase to shift the waveform by
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :return:
        waveform_time_domain: dict
            time-domain waveform polarisations
        waveform_frequency_domain: dict
            frequency-domain waveform polarisations
    """
    # Index shift
    waveform_time_domain = wrap_by_n_indices(int(index_shift), waveform_time_domain)
    # Phase shift
    waveform_frequency_domain = fourier_transform(
        waveform_time_domain, sampling_frequency
    )
    waveform_frequency_domain = {
        key: waveform_frequency_domain[key] * np.exp(-2j * phase_shift)
        for key in waveform_frequency_domain.keys()
    }
    waveform_time_domain = {
        key: bb.core.utils.infft(waveform_frequency_domain[key], sampling_frequency)
        for key in waveform_time_domain.keys()
    }
    return waveform_time_domain, waveform_frequency_domain


def calculate_overlaps_optimizable(new_params, *args):
    """
    Optimisable function for calculating overlaps.
    :param new_params: list
        list of new parameters to calculate the overlap with
                index_shift: int
                    number of indices to shift the waveform by
                phase_shift: float
                    phase to shift the waveform by
    :param args: list
        list of additional arguments
            waveform_time_domain: dict
                time-domain waveform polarisations
             comparison_waveform_frequency_domain: dict
                frequency-domain waveform polarisations for the comparison waveform
            frequency_array: array
                frequency array
            sampling_frequency: int
                frequency with which to 'sample' the waveform
            PSD: PowerSpectralDensity
                object to manage Power Spectral Densities
    :return:
        -overlap: float
            negative of the overlap between the two waveforms
    """
    # New guesses
    index_shift = new_params[0]
    phase_shift = new_params[1]
    # Extract arguments
    waveform_time_domain, comparison_waveform_frequency_domain, frequency_array, sampling_frequency, PSD = (
        args
    )
    # Apply shifts
    waveform_time_domain, waveform_frequency_domain = apply_shifts(
        waveform_time_domain, index_shift, phase_shift, sampling_frequency
    )
    return -overlap_function(
        a=waveform_frequency_domain,
        b=comparison_waveform_frequency_domain,
        frequency=frequency_array,
        psd=PSD,
    )


def maximise_overlap(
    waveform_time_domain,
    comparison_waveform_frequency_domain,
    sampling_frequency,
    frequency_array,
    PSD,
):
    """
    Funtion to determine the index and phase shifts that cause the SEOBNRe waveform to overlap
    maximally with a comparison waveform.
    :param waveform_time_domain: dict
        time-domain waveform polarisations
    :param comparison_waveform_frequency_domain: dict
        frequency-domain waveform polarisations of the comparison waveform
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :param frequency_array: array
        frequency array
    :param PSD: PowerSpectralDensity
        Object to manage the detector Power Spectral Densities
    :return:
        waveform_time_domain: dict
            time-domain waveform polarisations
        waveform_frequency_domain: dict
            frequency-domain waveform polarisations
        maximum_overlap: float
            value of maximum overlap, computed between the returned waveform and the comparison waveform
        index_shift: int
            the number of index shifts required to cause maximum overlap
        phase_shift: float
            the phase shift required to cause maximum overlap
    """
    comparison_waveform_time_domain = {
        key: bb.core.utils.infft(
            comparison_waveform_frequency_domain[key], sampling_frequency
        )
        for key in comparison_waveform_frequency_domain.keys()
    }
    comparison_max_index = np.argmax(
        np.abs(
            np.abs(comparison_waveform_time_domain["plus"])
            + np.abs(comparison_waveform_time_domain["cross"])
        )
    )
    # Apply a window
    waveform_time_domain = apply_tukey_window(waveform_time_domain, alpha=0.05)
    # Do the signal processing
    waveform_time_domain = process_signal(
        waveform_time_domain,
        comparison_length=len(comparison_waveform_time_domain["plus"]),
    )
    # Apply the initial time shift
    waveform_time_domain, initial_shift = wrap_at_maximum(
        waveform_time_domain, comparison_max_index
    )
    # Do the Fourier transform
    waveform_frequency_domain = fourier_transform(
        waveform_time_domain, sampling_frequency
    )
    # Now compute initial overlap
    maximum_overlap = overlap_function(
        comparison_waveform_frequency_domain,
        waveform_frequency_domain,
        frequency_array,
        PSD,
    )
    # Now we try to optimise this.
    index_limit = 0.005 * sampling_frequency
    phase_limit = np.pi / 2
    index_shift_guess = -index_limit
    phase_shift_guess = -phase_limit
    count = 0
    index_shift = index_shift_guess
    phase_shift = phase_shift_guess
    # First we hope to do this with the quick, dynamic method
    while np.round(maximum_overlap, 3) < 0.999:
        new_index_shift, new_phase_shift, new_overlap = fast_overlap_optimize(
            index_shift_guess,
            index_limit,
            phase_shift_guess,
            phase_limit,
            waveform_time_domain,
            comparison_waveform_frequency_domain,
            frequency_array,
            sampling_frequency,
            PSD,
        )
        index_shift_guess = random.random() * (
            index_shift_guess - 2 * index_shift_guess
        )
        phase_shift_guess = random.random() * (
            phase_shift_guess - 2 * phase_shift_guess
        )
        if new_overlap > maximum_overlap:
            maximum_overlap = new_overlap
            index_shift = new_index_shift
            phase_shift = new_phase_shift
        count += 1
        if count > 20:
            break
    # If that fails, we try the slow, rigorous method
    if np.round(maximum_overlap, 3) < 0.999:
        new_index_shift, new_phase_shift, new_overlap = slow_overlap_optimize(
            index_shift,
            index_limit,
            phase_shift,
            phase_limit,
            waveform_time_domain,
            comparison_waveform_frequency_domain,
            frequency_array,
            sampling_frequency,
            PSD,
            maximum_overlap=maximum_overlap,
        )
        if new_overlap > maximum_overlap:
            maximum_overlap = new_overlap
            index_shift = new_index_shift
            phase_shift = new_phase_shift
    # Generate the waveforms to return
    waveform_time_domain, waveform_frequency_domain = apply_shifts(
        waveform_time_domain, index_shift, phase_shift, sampling_frequency
    )
    return (
        waveform_time_domain,
        waveform_frequency_domain,
        maximum_overlap,
        index_shift,
        phase_shift,
    )


def slow_overlap_optimize(
    index_shift,
    index_limit,
    phase_shift,
    phase_limit,
    waveform_time_domain,
    comparison_waveform_frequency_domain,
    frequency_array,
    sampling_frequency,
    PSD,
    maximum_overlap=0,
):
    """
    Slow, rigorous, reliable method for determining the maximum overlap.
    :param index_shift: int
        suggested initial optimal index shift
    :param index_limit: int
        absolute value of the symmetric index shift limit
    :param phase_shift: float
        suggested initial optimal phase shift
    :param phase_limit:
        absolute value of the symmetric phase shift limit
    :param waveform_time_domain: dict
        time-domain waveform polarisations
    :param comparison_waveform_frequency_domain: dict
        frequency-domain waveform polarisations of the comparison waveform
    :param frequency_array: array
        frequency array
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :param PSD: PowerSpectralDensity
        object that manages the detector Power Spectral Densities
    :param maximum_overlap: float, 0
        current maximum overlap
    :return:
        index_shift: int
            index shift required for maximum overlap
        phase_shift: float
            phase shift required for maximum overlap
        maximum_overlap: float
            maximum overlap obtained
    """
    grid_size = 100
    # Flat time grids, no repeats
    phase_grid_init = np.linspace(-phase_limit, phase_limit, grid_size)
    time_grid_init = np.linspace(-index_limit, index_limit, grid_size)
    waveform_grid_time = [
        wrap_by_n_indices(int(n), waveform_time_domain) for n in time_grid_init
    ]
    waveform_grid_frequency = fourier_transform(waveform_grid_time, sampling_frequency)
    waveform_grid_shifted = [
        [
            {key: w[key] * np.exp(-2j * phase) for key in waveform_time_domain.keys()}
            for w in waveform_grid_frequency
        ]
        for phase in phase_grid_init
    ]
    overlap_grid = [
        [
            overlap_function(
                comparison_waveform_frequency_domain, w[t], frequency_array, PSD
            )
            for w in waveform_grid_shifted
        ]
        for t in range(grid_size)
    ]
    new_overlap = np.amax(overlap_grid)
    if new_overlap > maximum_overlap:
        maximum_overlap = new_overlap
        max_index = np.where(overlap_grid == maximum_overlap)
        max_index = max_index[0][0], max_index[1][0]
        index_shift = time_grid_init[max_index[0]]
        phase_shift = phase_grid_init[max_index[1]]
    return index_shift, phase_shift, maximum_overlap


def fast_overlap_optimize(
    index_shift_guess,
    index_limit,
    phase_shift_guess,
    phase_limit,
    waveform_time_domain,
    comparison_waveform_frequency_domain,
    frequency_array,
    sampling_frequency,
    PSD,
):
    """
    Fast, non-rigorous, but adequate method for estimating the maximum overlap.
    :param index_shift_guess: int
        suggested initial optimal index shift
    :param index_limit: int
        absolute value of the symmetric index shift limit
    :param phase_shift_guess: float
        suggested initial optimal phase shift
    :param phase_limit:
        absolute value of the symmetric phase shift limit
    :param waveform_time_domain: dict
        time-domain waveform polarisations
    :param comparison_waveform_frequency_domain: dict
        frequency-domain waveform polarisations of the comparison waveform
    :param frequency_array: array
        frequency array
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :param PSD: PowerSpectralDensity
        object that manages the detector Power Spectral Densities
    :return:
        index_shift: int
            index shift required for maximum overlap
        phase_shift: float
            phase shift required for maximum overlap
        maximum_overlap: float
            maximum overlap obtained
    """
    x0 = np.array([index_shift_guess, phase_shift_guess])
    bounds = [
        (index_shift_guess - index_limit, index_shift_guess + index_limit),
        (phase_shift_guess - phase_limit, phase_shift_guess + phase_limit),
    ]
    args = (
        waveform_time_domain,
        comparison_waveform_frequency_domain,
        frequency_array,
        sampling_frequency,
        PSD,
    )
    res = minimize(
        calculate_overlaps_optimizable, bounds=bounds, x0=x0, args=args, tol=1e-100
    )
    index_shift, phase_shift = res.x[0], res.x[1]
    maximum_overlap = -res.fun
    return index_shift, phase_shift, maximum_overlap

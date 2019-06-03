"""

A file for generating waveforms and related objects for the pySEOBNRe package.

"""

import subprocess
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import bilby as bb


def read_in_seobnre(filename):
    """
    Read in an SEOBNRe output file of a waveform strain time series.
    :param filename: str
        the name of the SEOBNRe output file
    :return:
        t: array
            time array
        waveform: dict
            time-domain waveform polarisations
    """
    t = []
    hp = []
    hc = []
    for line in open(filename, "r"):
        if "#" not in line:
            line_split = line.split()
            t.append(float(line_split[0]))
            hp.append(float(line_split[1]))
            hc.append(float(line_split[2]))
    return t, dict(plus=hp, cross=hc)


def seobnre_bbh_with_spin_and_eccentricity(
    parameters, sampling_frequency, minimum_frequency, maximum_frequency=0, plot=False
):
    """
    Return the time array and waveform polarisations simulated by the SEOBNRe c-code.
    :param parameters: dict
        dictionary of waveform parameters
    :param sampling_frequency: int
        frequency with which to 'sample' the waveform
    :param minimum_frequency: int
        minimum frequency to contain in the signal
    :param maximum_frequency: int, 0
        maximum frequency to contain in the signal. if 0, create as much as possible
    :param plot: Boolean, False
        if True, plot the generated waveform
    :return:
        t: array
            time array
        seobnre: dict
            time-domain waveform polarisations
    """
    make = ["make"]
    c_code = "/Users/irom0001/eccentricity/SEOBNRE/package/c_code/"
    phiRef = parameters["phase"]
    m1 = parameters["mass_1"]
    m2 = parameters["mass_2"]
    f_min = minimum_frequency
    f_max = maximum_frequency
    f_samp = sampling_frequency
    e0 = parameters["eccentricity"]
    distance = parameters["luminosity_distance"]
    inclination = parameters["theta_jn"]
    s1z = parameters["chi_1"]
    s1y = 0
    s1x = 0
    s2z = parameters["chi_2"]
    s2y = 0
    s2x = 0
    outfile_name = "simulation_" + "{}_{}_{}_{}_{}_{}_{}.dat".format(
        m1, m2, distance, e0, inclination, s1z, s2z
    )
    # Generate the SEOBNRe time-domain waveform
    execute = (
        "./SEOBNRE "
        "--phiRef {} --m1 {} --m2 {} --e0 {} --distance {} --inclination {} "
        "--spin1x {} --spin1y {} --spin1z {} --spin2x {} --spin2y {} --spin2z {} "
        "--f-min {} --f-max {} --sample-rate {} "
        "--outname {}".format(
            phiRef,
            m1,
            m2,
            e0,
            distance,
            inclination,
            s1x,
            s1y,
            s1z,
            s2x,
            s2y,
            s2z,
            f_min,
            f_max,
            f_samp,
            outfile_name,
        ).split()
    )
    # Make
    subprocess.Popen(
        make, cwd=c_code, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).wait()
    # Execute
    subprocess.Popen(
        execute, cwd=c_code, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).wait()
    # Read in the time-domain waveform
    t, seobnre = read_in_seobnre(c_code + outfile_name)
    # Plot if requested
    if plot:
        fig = plt.figure()
        plt.plot(t, seobnre["plus"], label="SEOBNRe+")
        plt.plot(t, seobnre["cross"], label="SEOBNRex")
        plt.xlabel("Time (s)")
        plt.ylabel("Strain")
        plt.legend()
        plt.show()
    # Clear up the file
    subprocess.Popen(
        ["rm", outfile_name], cwd=c_code, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return t, seobnre


def get_IMRPhenomD_comparison_waveform_generator(
    minimum_frequency, sampling_frequency, duration
):
    """
    Provide the waveform generator object for the comparison waveform, IMRPhenomD.
    :param minimum_frequency: int
        minimum frequency to contain in the waveform
    :param sampling_frequency: int
        frequency with which to 'sample' the signal
    :param duration: int
        time duration of the signal
    :return:
        waveform_generator: WaveformGenerator
            the waveform generator object for the comparison waveform
    """
    waveform_approximant = "IMRPhenomD"
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant,
        reference_frequency=minimum_frequency,
        minimum_frequency=minimum_frequency,
    )
    waveform_generator = bb.gw.WaveformGenerator(
        duration=duration,
        sampling_frequency=sampling_frequency,
        frequency_domain_source_model=bb.gw.source.lal_binary_black_hole,
        parameter_conversion=bb.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_arguments,
    )
    return waveform_generator

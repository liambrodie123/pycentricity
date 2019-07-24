"""

A file for reweighting non-SEOBNRe results using an SEOBNRe waveform for the pySEOBNRe package.

"""

import bilby as bb
import numpy as np
import python_code.overlap as ovlp
import python_code.waveform as wf
import numpy.random as random

np.random.seed(3)
minimum_eccentricity = 1e-4


def log_likelihood_interferometer(
    waveform_polarizations, interferometer, parameters, duration
):
    """
    Return the log likelihood at a single interferometer.
    :param waveform_polarizations: dict
        frequency-domain waveform polarisations
    :param interferometer: Interferometer
        Interferometer object
    :param parameters: dict
        waveform parameters
    :param duration: int
        time duration of the signal
    :return:
        log_l: float
            log likelihood evaluation
    """
    signal_ifo = interferometer.get_detector_response(
        waveform_polarizations, parameters
    )
    log_l = (
        -2.0
        / duration
        * np.vdot(
            interferometer.frequency_domain_strain - signal_ifo,
            (interferometer.frequency_domain_strain - signal_ifo)
            / interferometer.power_spectral_density_array,
        )
    )
    return log_l.real


def noise_log_likelihood(interferometers, duration):
    """
    Return the likelihood that the data is caused by noise
    :param interferometers: InterferometerList
        list of interferometers involved in the detection
    :param duration: int
        time duration of the signal
    :return:
        log_l: float
        log likelihood evaluation
    """
    log_l = 0
    for interferometer in interferometers:
        log_l -= (
            2.0
            / duration
            * np.sum(
                abs(interferometer.frequency_domain_strain) ** 2
                / interferometer.power_spectral_density_array
            )
        )
    return log_l.real


def log_likelihood(waveform_polarizations, interferometers, parameters, duration):
    """
    Return the log likelihood of all interferometers involved in the detection
    :param waveform_polarizations: dict
        frequency-domain waveform polarisations
    :param interferometers: InterferometerList
        list of interferometers involved in the detection
    :param parameters: dict
        waveform parameters
    :param duration: int
        time duration of the signal
    :return:
        log_l: float
        log likelihood evaluation
    """
    log_l = 0
    for interferometer in interferometers:
        log_l += log_likelihood_interferometer(
            waveform_polarizations, interferometer, parameters, duration
        )
    return log_l.real


def log_likelihood_ratio(waveform_polarizations, interferometers, parameters, duration):
    """
    Return the ratio of the likelihood of the signal being real to the likelihood
    of the signal being noise
    :param waveform_polarizations: dict
        frequency-domain waveform polarisations
    :param interferometers: InterferometerList
        list of interferometers involved in the relationship
    :param parameters: dict
        waveform variables
    :param duration: int
        time duration of the signal
    :return:
        log_likelihood_ratio: float
            log likelihood evaluation ratio
    """
    return log_likelihood(
        waveform_polarizations, interferometers, parameters, duration
    ) - noise_log_likelihood(interferometers, duration)


def pick_weighted_random_eccentricity(cumulative_density_grid, eccentricity_grid):
    """
    Return a random eccentricity, weighted by a cumulative density function.
    :param cumulative_density_grid: array
        1D grid of cumulative densities
    :param eccentricity_grid: array
        1D grid of eccentricities
    :return:
        random_eccentricity: float
            the weighted random eccentricity chosen
    """
    # First select a bin, weighted at random
    start = cumulative_density_grid[0]
    end = cumulative_density_grid[-1]
    random_value = random.random_sample() * (end - start) + start
    for i, cd in enumerate(cumulative_density_grid):
        if cd >= random_value:
            # Now select an eccentricity at random from within the selected bin
            lower_bound = 0
            upper_bound = eccentricity_grid[i]
            if i > 0:
                lower_bound = eccentricity_grid[i - 1]
            random_eccentricity = (
                random.random_sample() * (upper_bound - lower_bound) + lower_bound
            )
            return random_eccentricity


def cumulative_density_function(log_likelihood_grid):
    """
    Return a cumulative density function over a grid of eccentricities.
    :param log_likelihood_grid: array
        1D grid of log likelihoods
    :return:
    cumulative_density: array
        1D grid of cumulative densities
    """
    # IRS - (temporary?) edits to deal with extremely high values of log likelihood
    minimum_log_likelihood = np.min(log_likelihood_grid)
    # Ratio of likelihood to minimum log likelihood
    log_likelihood_grid = log_likelihood_grid - minimum_log_likelihood
    likelihood_grid = np.exp(log_likelihood_grid)
    cumulative_density = np.cumsum(likelihood_grid)
    cumulative_density_normalised = cumulative_density / cumulative_density[-1]
    return cumulative_density_normalised


def new_weight(
    log_L,
    parameters,
    comparison_waveform_frequency_domain,
    interferometers,
    duration,
    sampling_frequency,
    maximum_frequency,
    label,
):
    """
    Compute the new weight for a point, weighted by the eccentricity-marginalised likelihood.
    :param log_L: float
        the original likelihood of the point
    :param parameters: dict
        the parameters that define the sample
    :param comparison_waveform_frequency_domain: dict
        frequency-domain waveform polarisations of the waveform used for the original analysis
    :param interferometers: InterferometerList
        list of interferometers used in the detection
    :param duration: int
        time duration of the signal
    :param sampling_frequency: int
        the frequency with which to 'sample' the waveform
    :param minimum_frequency: int
        the minimum frequency at which the data is analysed
    :param maximum_frequency: int
        the maximum frequency at which the data is analysed
    :param label: str
        identifier for results
    :return:
        e: float
            the new eccentricity sample
        average_log_likelihood: float
            the eccentricity-marginalised new likelihood
        log_weight: float
            the log weight of the sample
    """
    # First calculate a grid of likelihoods.
    grid_size = 20
    eccentricity_grid = np.logspace(
        np.log10(minimum_eccentricity), np.log10(0.2), grid_size
    )
    # Recalculate the log likelihood of the original sample
    recalculated_log_likelihood = log_likelihood_ratio(
        comparison_waveform_frequency_domain, interferometers, parameters, duration
    )
    # Print a warning if this is much different to the likelihood stored in the results
    if abs(recalculated_log_likelihood - log_L) / log_L > 0.1:
        percentage = abs(recalculated_log_likelihood - log_L) / log_L * 100
        print(
            "WARNING :: recalculated log likelihood differs from original by {}%".format(
                percentage
            )
        )
        print('original log L: ' + str(log_L))
        print('recalculated log L: ' + str(recalculated_log_likelihood))
    log_likelihood_grid = []
    intermediate_outfile = open(label + "_eccentricity_result.txt", "w")
    intermediate_outfile.write("sample parameters:\n")
    for key in parameters.keys():
        intermediate_outfile.write(key + ":\t" + str(parameters[key]) + "\n")
    intermediate_outfile.write("\n-------------------------\n")
    intermediate_outfile.write("e\t\tlog_L\t\tmaximised_overlap\n")
    # Prepare for the possibility that we have to disregard this sample
    disregard = False
    for e in eccentricity_grid:
        parameters.update({"eccentricity": e})
        # Need to have a set minimum frequency, since this is also the reference frequency
        t, seobnre_waveform_time_domain = wf.seobnre_bbh_with_spin_and_eccentricity(
            parameters=parameters,
            sampling_frequency=sampling_frequency,
            minimum_frequency=10,
            maximum_frequency=maximum_frequency + 1000,
        )
        if t is None:
            print('No waveform generated; disregard sample ' + label)
            intermediate_outfile.write(
                str(e)
                + "\t\t"
                + str(None)
                + "\t\t"
                + str(None)
                + "\n"
            )
            disregard = True
        else:
            seobnre_wf_td, seobnre_wf_fd, max_overlap, index_shift, phase_shift = ovlp.maximise_overlap(
                seobnre_waveform_time_domain,
                comparison_waveform_frequency_domain,
                sampling_frequency,
                interferometers[0].frequency_array,
                interferometers[0].power_spectral_density,
            )
            seobnre_wf_fd = ovlp.zero_pad_frequency_domain_signal(
                seobnre_wf_fd, interferometers
            )
            eccentric_log_L = log_likelihood_ratio(
                seobnre_wf_fd, interferometers, parameters, duration
            )
            log_likelihood_grid.append(eccentric_log_L)
            intermediate_outfile.write(
                str(e)
                + "\t\t"
                + str(eccentric_log_L)
                + "\t\t"
                + str(max_overlap)
                + "\n"
            )
    intermediate_outfile.close()
    if not disregard:
        cumulative_density_grid = cumulative_density_function(log_likelihood_grid)
        # We want to pick a weighted random point from within the CDF
        e = pick_weighted_random_eccentricity(cumulative_density_grid, eccentricity_grid)
        # Also return eccentricity-marginalised log-likelihood
        average_log_likelihood = np.mean(log_likelihood_grid)
        # The weight is the ratio of this to the log likelihood
        log_weight = average_log_likelihood - recalculated_log_likelihood
        return e, average_log_likelihood, log_weight
    else:
        return None, None, None


def reweight_by_eccentricity(
    samples,
    log_likelihood,
    sampling_frequency,
    comparison_waveform_generator,
    interferometers,
    duration,
    output_folder=".",
    maximum_frequency=0,
    label="",
):
    """
    Function to return a dictionary containing the eccentricity-marginalised log likelihood,
    the weight of the sample and the new eccentricity.
    :param samples: dict
        dictionary of all posterior samples from the original analysis
    :param log_likelihood: list
        list of log-likelihoods from the original analysis
    :param sampling_frequency: int
        the frequency at which the data is sampled
    :param minimum_frequency: int
        the minimum frequency at which the data is sampled
    :param comparison_waveform_generator: WaveformGenerator
        the waveform generator for the waveform used in the original analysis
    :param interferometers: InterferometerList
        list of interferometers used in the detection
    :param duration: int
        time duration of the analysis
    :param output_folder: str
        location to send the output
    :param maximum_frequency: int
        the maximum frequency of the analysis
    :param label: str
        identifier for results
    :return:
        output: dict
            dictionary of output from the reweighting procedure
    """
    number_of_samples = len(log_likelihood)
    converted_samples, added_keys = bb.gw.conversion.convert_to_lal_binary_black_hole_parameters(
        samples
    )
    converted_samples = {
        key: converted_samples[key].values.tolist()
        for key in converted_samples.keys()
        if type(converted_samples[key]) is not float
    }
    # List of parameters - this will be what we use to compute the highest overlap
    parameter_list = [
        dict(
            mass_1=converted_samples["mass_1"][i][0],
            mass_2=converted_samples["mass_2"][i][0],
            luminosity_distance=converted_samples["luminosity_distance"][i][0],
            geocent_time=converted_samples["geocent_time"][i][0],
            ra=converted_samples["ra"][i][0],
            dec=converted_samples["dec"][i][0],
            theta_jn=converted_samples["theta_jn"][i][0],
            psi=converted_samples["psi"][i][0],
            phase=converted_samples["phase"][i][0],
            chi_1=converted_samples["chi_1"][i][0],
            chi_2=converted_samples["chi_2"][i][0],
            eccentricity=minimum_eccentricity,
        )
        for i in range(number_of_samples)
    ]
    # Generate the IMRPhenomD waveform for each sample
    comparison_waveform_strain_list = [
        comparison_waveform_generator.frequency_domain_strain(parameters)
        for parameters in parameter_list
    ]
    output = {key: [] for key in ["eccentricity", "new_log_L", "log_weight"]}
    # Write the output file along the way
    outfile = open(output_folder + "/" + label + "_master_output_store.txt", "w")
    outfile.write("i\t\te\t\tnew_log_L\t\tlog_w\n")
    for i, log_L in enumerate(log_likelihood):
        # If the spins are too large, the sample may fail to generate eccentric waveforms,
        # so we impose a moderate-spin prior here
        if any([parameter_list[i]['chi_1'] > 0.6, parameter_list[i]['chi_2'] > 0.6]):
            print(
                    'omitting sample; chi_1 = ' + str(parameter_list[i]['chi_1'])
                    + ', chi_2 = ' + str(parameter_list[i]['chi_2'])
            )
            output["eccentricity"].append(None)
            output["new_log_L"].append(None)
            output["log_weight"].append(None)
            outfile.write(
                str(i)
                + "\t\t"
                + str(None)
                + "\t\t"
                + str(None)
                + "\t\t"
                + str(None)
                + "\n"
            )
        else:
            eccentricity, new_log_L, log_weight = new_weight(
                log_L,
                parameter_list[i],
                comparison_waveform_strain_list[i],
                interferometers,
                duration,
                sampling_frequency,
                maximum_frequency,
                output_folder + "/" + label + "_" + str(i),
            )
            outfile.write(
                str(i)
                + "\t\t"
                + str(eccentricity)
                + "\t\t"
                + str(new_log_L)
                + "\t\t"
                + str(log_weight)
                + "\n"
            )
            output["eccentricity"].append(eccentricity)
            output["new_log_L"].append(new_log_L)
            output["log_weight"].append(log_weight)
        print(
            "new weight calculation {}% complete".format(
                np.round(i / number_of_samples * 100, 2)
            )
        )
    outfile.close()
    return output

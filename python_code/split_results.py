"""

A file for splitting a results file into the required subset for the pySEOBNRe package.

"""

import bilby as bb
import numpy as np
import json


def result_subset(start_index, end_index, result):
    """
    Return a certain subset of a set of results.
    :param start_index: int
        the index at which to begin the subset
    :param end_index: int
        the index at which to end the subset
    :param result: Result
        the result object
    :return:
        results_subset_dictionary: dict
            dictionary containing information about the result subset
    """
    samples = result.posterior
    log_likelihood = result.posterior.log_likelihood
    subset_samples = {
        key: samples[key][start_index:end_index] for key in samples.keys()
    }
    subset_log_likelihood = log_likelihood[start_index:end_index]
    results_subset_dictionary = dict(
        samples=subset_samples, log_likelihoods=subset_log_likelihood
    )
    return results_subset_dictionary


def write_subset_to_file(start_index, end_index, result, output_file_path):
    """
    Save a certain subset to a file.
    :param start_index: int
        the index at which to begin the subset
    :param end_index: int
        the index at which to end the subset
    :param result: Result
        the result object
    """
    result_subset_dictionary = result_subset(start_index, end_index, result)
    output_file_path += "/result_subset_s{}_e{}.json".format(start_index, end_index)
    with open(output_file_path, "w") as f:
        json.dump(result_subset_dictionary, f)


def split_results_into_subsets(number_per_file, result_file):
    # Work out where to store the output
    output_file_path_list = result_file.split("/")
    output_file_path = ""
    for string in output_file_path_list[0:-1]:
        output_file_path += string
    # Get the result object
    result = bb.result.read_in_result(result_file)
    total_number_of_samples = len(result.posterior.log_likelihood)
    print(total_number_of_samples)
    start_indices = np.arange(0, total_number_of_samples, number_per_file)
    print(start_indices)
    print(len(start_indices))
    end_indices = np.arange(number_per_file, total_number_of_samples, number_per_file)
    end_indices = np.concatenate((end_indices, [total_number_of_samples]))
    print(end_indices)
    print(len(end_indices))


split_results_into_subsets(
    50,
    "/Users/irom0001/eccentricity/SEOBNRE/package/bilby_results/GW150914/dynesty_GW150914_IMRPhenomD_BWpsd_for_eccentricity_dynesty_fix_combined_result.json",
)

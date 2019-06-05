"""

A file for splitting a results file into the required subset for the pySEOBNRe package.

"""

import python_code.utils as utils
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
    keys = utils.search_keys
    log_likelihood = result.posterior.log_likelihood
    subset_samples = {key: samples[key][start_index:end_index].tolist() for key in keys}
    subset_log_likelihood = log_likelihood[start_index:end_index].tolist()
    results_subset_dictionary = dict(
        samples=subset_samples, log_likelihoods=subset_log_likelihood
    )
    return results_subset_dictionary


def write_subset_to_file(start_index, end_index, result, index, output_file_path):
    """
    Save a certain subset to a file.
    :param start_index: int
        the index at which to begin the subset
    :param end_index: int
        the index at which to end the subset
    :param result: Result
        the result object
    :param output_file_path: str
        the path to the location in which to create output files
    """
    result_subset_dictionary = result_subset(start_index, end_index, result)
    output_file_path += "/result_{}.json".format(index)
    with open(output_file_path, "w") as f:
        json.dump(result_subset_dictionary, f)


def split_results_into_subsets(number_per_file, result_file):
    """
    Split a set of results into numerous subsets for easier computation
    :param number_per_file: int
        number of results to store per file
    :param result_file: str
        file path of results file
    """
    # Work out where to store the output
    output_file_path_list = result_file.split("/")
    output_file_path = ""
    for string in output_file_path_list[0:-1]:
        output_file_path += string + "/"
    output_file_path += "subsets/"
    bb.core.utils.check_directory_exists_and_if_not_mkdir(output_file_path)
    # Get the result object
    result = bb.result.read_in_result(result_file)
    total_number_of_samples = len(result.posterior.log_likelihood)
    start_indices = np.arange(0, total_number_of_samples, number_per_file)
    for i, start_index in enumerate(start_indices):
        end_index = min(start_index + number_per_file, total_number_of_samples)
        write_subset_to_file(start_index, end_index, result, i, output_file_path)
    print("results split into {} subsets".format(i))

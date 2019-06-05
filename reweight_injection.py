import bilby as bb
import pandas as pd

import python_code.reweight as rwt
import python_code.waveform as wf

import json
import argparse


# Set up the argument parser
parser = argparse.ArgumentParser("Produce weights for one subset of a results file")
parser.add_argument("-s", "--sub-result", help="Path of the sub-result file to use")
args = parser.parse_args()

# Access the samples
json_data = json.load(open(args.sub_result))
samples = json_data["samples"]
samples = {key: pd.DataFrame(samples[key]) for key in samples.keys()}
log_likelihoods = pd.DataFrame(json_data["log_likelihoods"])

# Set up the basic properties of the runs
maximum_frequency = 1024
sampling_frequency = 4096
post_trigger_duration = 2
deltaT = 0.2

# Read event-specific properties from the utils file
minimum_frequency = 20
duration = 8
detectors = ['H1', 'L1']

# Generate the comparison waveform generator
waveform_generator = wf.get_IMRPhenomD_comparison_waveform_generator(
    minimum_frequency, sampling_frequency, duration
)
# Frequency array
frequency_array = waveform_generator.frequency_array

# Set up the interferometers with event data
interferometers = bb.gw.detector.InterferometerList(detectors)
interferometers.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=-6,
)
for ifo in interferometers:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency

# Output to the folder with all of the result subsets
folder_list = args.sub_result.split("/")
folder = ""
for string in folder_list[0:-1]:
    folder += string + "/"
folder += "weights/"
bb.core.utils.check_directory_exists_and_if_not_mkdir(folder)
label = folder_list[-1].split(".")[0]
output = rwt.reweight_by_eccentricity(
    samples,
    log_likelihoods,
    sampling_frequency,
    minimum_frequency,
    waveform_generator,
    interferometers,
    duration,
    folder,
    maximum_frequency,
    label=label,
)
print("Results weighted for file " + args.sub_result)
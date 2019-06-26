#!/home/isobel.romero-shaw/anaconda3/bin/python3.6
import bilby as bb
import pandas as pd

import python_code.utils as utils
import python_code.reweight as rwt
import python_code.waveform as wf
import python_code.overlap as ovlp

import numpy as np

import json
import argparse


# Set up the argument parser
parser = argparse.ArgumentParser("Produce zero-eccentricity weights for one subset of a results file")
parser.add_argument("-e", "--event", help="GW-name of the event")
parser.add_argument("-s", "--sub-result", help="Path of the sub-result file to use")
args = parser.parse_args()

# Access the samples
json_data = json.load(open(args.sub_result))
samples = json_data["samples"]
samples = {key: pd.DataFrame(samples[key]) for key in samples.keys()}
log_likelihoods = json_data["log_likelihoods"]

# Set up the basic properties of the runs
maximum_frequency = 1024
post_trigger_duration = 2
deltaT = 0.2
wf_minimum_frequency = 20

# Read event-specific properties from the utils file
sampling_frequency = utils.sampling_frequency[args.event]
ifo_minimum_frequency = utils.minimum_frequency[args.event]
duration = utils.event_duration[args.event]
detectors = utils.event_detectors[args.event]
trigger_time = utils.trigger_time[args.event]

# Generate the comparison waveform generator
waveform_generator = wf.get_IMRPhenomD_comparison_waveform_generator(
    wf_minimum_frequency, sampling_frequency, duration
)
# Frequency array
frequency_array = waveform_generator.frequency_array

# Set up the interferometers with event data
interferometers = bb.gw.detector.InterferometerList(detectors)
start = trigger_time + post_trigger_duration - duration
end = start + duration
channel_dict = {
    key: key + ":" + utils.event_channels[args.event][key] for key in detectors
}
for ifo in interferometers:
    ifo.set_strain_data_from_csv(
        '/home/isobel.romero-shaw/public_html/PYCENTRICITY/pycentricity/submissions/'
        + args.event + '/event_data/' + ifo.name + '_time_domain_strain_data.csv'
    )
    ifo.power_spectral_density = bb.gw.detector.PowerSpectralDensity.from_power_spectral_density_file(
        psd_file=utils.event_psd_file_path[args.event][ifo.name]
    )
    ifo.minimum_frequency = ifo_minimum_frequency
    ifo.maximum_frequency = maximum_frequency

# Output to the folder with all of the result subsets
folder_list = args.sub_result.split("/")
folder = ""
for string in folder_list[0:-2]:
    folder += string + "/"
folder += "zero_eccentricity_weights/"
bb.core.utils.check_directory_exists_and_if_not_mkdir(folder)
label = folder_list[-1].split(".")[0]

number_of_samples = len(log_likelihoods)
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
        eccentricity=0.0,
    )
    for i in range(number_of_samples)
]
# Generate the IMRPhenomD waveform for each sample
comparison_waveform_strain_list = [
    waveform_generator.frequency_domain_strain(parameters)
    for parameters in parameter_list
]
output = {key: [] for key in ["eccentricity", "new_log_L", "log_weight"]}
# Write the output file along the way
outfile = open(folder + "/" + label + "_master_output_store.txt", "w")
outfile.write("i\t\te\t\tnew_log_L\t\tlog_w\n")
for i, log_L in enumerate(log_likelihoods):
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
        # Prepare for the possibility that we have to disregard this sample
        disregard = False
        # Need to have a set minimum frequency, since this is also the reference frequency
        t, seobnre_waveform_time_domain = wf.seobnre_bbh_with_spin_and_eccentricity(
            parameters=parameter_list[i],
            sampling_frequency=sampling_frequency,
            minimum_frequency=10,
            maximum_frequency=maximum_frequency + 1000,
        )
        if t is None:
            print('No waveform generated; disregard sample ' + label)
            disregard = True
        else:
            seobnre_wf_td, seobnre_wf_fd, max_overlap, index_shift, phase_shift = ovlp.maximise_overlap(
                seobnre_waveform_time_domain,
                comparison_waveform_strain_list[i],
                sampling_frequency,
                interferometers[0].frequency_array,
                interferometers[0].power_spectral_density,
            )
            seobnre_wf_fd = ovlp.zero_pad_frequency_domain_signal(
                seobnre_wf_fd, interferometers
            )
            new_log_L = rwt.log_likelihood_ratio(
                seobnre_wf_fd, interferometers, parameter_list[i], duration
            )
        if not disregard:
            # We want to pick a weighted random point from within the CDF
            eccentricity = 0
            # The weight is the ratio of this to the log likelihood
            recalculated_log_likelihood = rwt.log_likelihood_ratio(
                comparison_waveform_strain_list[i], interferometers, parameter_list[i], duration
            )
            log_weight = new_log_L - recalculated_log_likelihood
        else:
            eccentricity = None
            new_log_L = None
            log_weight = None

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

print("Results weighted for file " + args.sub_result + " for event " + args.event)
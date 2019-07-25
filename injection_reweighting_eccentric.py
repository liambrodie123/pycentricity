#!/home/isobel.romero-shaw/anaconda3/bin/python3.6
import bilby as bb
import python_code.waveform as wf
import python_code.overlap as ovlp
import python_code.reweight as rwt
import numpy as np
import json
import pandas as pd
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser("")
parser.add_argument("-i", "--index", help="")
args = parser.parse_args()

np.random.seed(5432)
resultdir = '/home/isobel.romero-shaw/public_html/PYCENTRICITY/pycentricity/injection_recovery_output/eccentric/subsets/'
weightsdir = '/home/isobel.romero-shaw/public_html/PYCENTRICITY/pycentricity/injection_recovery_output/eccentric/weights/'

# Access the samples
sub_result = resultdir + 'result_' + str(args.index) + '.json'
json_data = json.load(open(args.sub_result))
samples = json_data["samples"]
samples = {key: pd.DataFrame(samples[key]) for key in samples.keys()}
log_likelihoods = json_data["log_likelihoods"]

# Make the weights folder
bb.core.utils.check_directory_exists_and_if_not_mkdir(weightsdir)
label = 'eccentric_reweighting'

# injection parameters
injection_parameters = dict(
    mass_1=35.0,
    mass_2=30.0,
    eccentricity=0.1,
    luminosity_distance=440.0,
    theta_jn=0.4,
    psi=0.1,
    phase=1.2,
    geocent_time=0.0,
    ra=3.7,
    dec=5.73,
    chi_1=0.0,
    chi_2=0.0,
)

# Frequency settings
minimum_frequency = 20
maximum_frequency = 2046
sampling_frequency = 4096
duration = 8
deltaT = 0.2

# Comparison waveform
comparison_waveform_generator = wf.get_IMRPhenomD_comparison_waveform_generator(
    minimum_frequency, sampling_frequency, duration
)
comparison_waveform_frequency_domain = comparison_waveform_generator.frequency_domain_strain(
    injection_parameters
)
# Interferometers
interferometers = bb.gw.detector.InterferometerList(["H1", "L1"])
interferometers.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - duration + 2,
)
for ifo in interferometers:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency

# SEBONRe waveform to inject
t, seobnre_waveform_time_domain = wf.seobnre_bbh_with_spin_and_eccentricity(
    parameters=injection_parameters,
    sampling_frequency=sampling_frequency,
    minimum_frequency=minimum_frequency - 10,
    maximum_frequency=maximum_frequency + 1000,
)
seobnre_wf_td, seobnre_wf_fd, max_overlap, index_shift, phase_shift = ovlp.maximise_overlap(
    seobnre_waveform_time_domain,
    comparison_waveform_frequency_domain,
    sampling_frequency,
    interferometers[0].frequency_array,
    interferometers[0].power_spectral_density,
)
print('maximum overlap: ' + str(max_overlap))

seobnre_wf_fd = ovlp.zero_pad_frequency_domain_signal(seobnre_wf_fd, interferometers)
# Inject the signal
interferometers.inject_signal(
    parameters=injection_parameters, injection_polarizations=seobnre_wf_fd
)

# Output to the folder with all of the result subsets
folder_list = args.sub_result.split("/")
folder = ""
for string in folder_list[0:-2]:
    folder += string + "/"
folder += "weights/"
bb.core.utils.check_directory_exists_and_if_not_mkdir(folder)
label = folder_list[-1].split(".")[0]
output = rwt.reweight_by_eccentricity(
        samples,
        log_likelihoods,
        sampling_frequency,
        comparison_waveform_generator,
        interferometers,
        duration,
        folder,
        maximum_frequency,
        label=label,
)
print("Results weighted for file " + args.sub_result)






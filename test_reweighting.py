import bilby as bb
import python_code.reweight as rwt
import python_code.waveform as wf
import python_code.utils as utils
import gwpy


date = '150914'
event = 'GW' + date
folder = "/home/isobel.romero-shaw/public_html/IMRPD_for_reweighting/" + date + "/results_dynesty_fix/result/"
result = bb.result.read_in_result(
    folder
    + "dynesty_" + event + "_IMRPhenomD_BWpsd_for_eccentricity_dynesty_fix_combined_result.json"
)
samples = result.posterior
log_likelihood = result.posterior.log_likelihood
number_of_samples = len(log_likelihood)
number_to_test = 3
start_index = 3
test_samples = {
    key: samples[key][0:number_to_test] for key in samples.keys()
}
test_log_likelihood = log_likelihood[0:number_to_test]
minimum_frequency = utils.minimum_frequency[event]
reference_frequency = 20
maximum_frequency = 1024
sampling_frequency = 4096
duration = utils.event_duration[event]
post_trigger_duration = 2
n_nodes = 5
deltaT = 0.2

detectors = utils.event_detectors[event]
trigger_time = utils.trigger_time[event]


waveform_generator = wf.get_IMRPhenomD_comparison_waveform_generator(
    minimum_frequency, sampling_frequency, duration
)
frequency_array = waveform_generator.frequency_array

ifos = bb.gw.detector.InterferometerList(detectors)
start = trigger_time + post_trigger_duration - duration
end = start + duration
channel_dict = {
    key: key + ':' + utils.event_channels[event][key] for key in detectors
}
for ifo in ifos:
    data = gwpy.timeseries.TimeSeries.get(
        channel_dict[ifo.name], start, end, verbose=False
    )
    data = data.resample(sampling_frequency)
    ifo.strain_data.set_from_gwpy_timeseries(data)
    ifo.power_spectral_density = bb.gw.detector.PowerSpectralDensity.from_power_spectral_density_file(
        psd_file=utils.event_psd_file_path[event][ifo.name]
    )
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency

output = rwt.reweight_by_eccentricity(
    test_samples,
    test_log_likelihood,
    sampling_frequency,
    minimum_frequency,
    waveform_generator,
    ifos,
    duration,
    folder,
    maximum_frequency
)
print(output)

"""
    The user must run this script prior to submitting a job in order to point gwpy.TimeSeries.read to data
    and avoid downloading the same data for every job
"""

import bilby as bb
import pandas as pd

import python_code.utils as utils
import python_code.reweight as rwt
import python_code.waveform as wf

import json
import gwpy
import argparse


# Set up the argument parser
parser = argparse.ArgumentParser("Download the data for an event")
parser.add_argument("-e", "--event", help="GW-name of the event")
args = parser.parse_args()
# Outdir
outdir = "submissions/" + args.event + "/event_data/"
# Get the data for the right segment
duration = utils.event_duration[args.event]
detectors = utils.event_detectors[args.event]
trigger_time = utils.trigger_time[args.event]

post_trigger_duration = 2
sampling_frequency = 4096

start = trigger_time + post_trigger_duration - duration
end = start + duration
# Set up the interferometers with event data
interferometers = bb.gw.detector.InterferometerList(detectors)
channel_dict = {
    key: key + ":" + utils.event_channels[args.event][key] for key in detectors
}
bb.core.utils.check_directory_exists_and_if_not_mkdir(outdir)
for ifo in interferometers:
    data = gwpy.timeseries.TimeSeries.get(
        channel_dict[ifo.name], start, end, verbose=False
    )
    data = data.resample(sampling_frequency)
    ifo.strain_data.set_from_gwpy_timeseries(data)
    with open(outdir + '/' + ifo.name + '_time_domain_strain_data.csv', 'w') as f:
        for i, t in enumerate(ifo.time_array):
            f.write(str(t) + ',' + str(ifo.time_domain_strain[i]) + '\n')

interferometers.save_data(outdir=outdir, label=args.event)
interferometers.plot_data()


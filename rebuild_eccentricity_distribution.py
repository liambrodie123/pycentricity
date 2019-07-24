import python_code.reweight as rwt

import glob
import argparse

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams


fontparams = {'mathtext.fontset': 'stix',
             'font.family': 'serif',
             'font.serif': "Times New Roman",
             'mathtext.rm': "Times New Roman",
             'mathtext.it': "Times New Roman:italic",
             'mathtext.sf': 'Times New Roman',
             'mathtext.tt': 'Times New Roman'}
rcParams.update(fontparams)

# Set up the argument parser
parser = argparse.ArgumentParser("Rebuild the eccentricity distribution based on weights given in output files")
parser.add_argument("-r", "--results-folder", help="Location of the results files")
args = parser.parse_args()

result_files = list(glob.glob(args.results_folder + '/result_*_*_eccentricity_result.txt'))
eccentricity_distribution = []

for single_file in result_files:
    disregard = False
    with open(single_file, 'r') as s_f:
        read_from_here = False
        log_L = []
        eccentricities = []
        for line in s_f:
            if read_from_here:
                if 'None' not in line:
                    split_line = line.split('\t\t')
                    eccentricities.append(float(split_line[0]))
                    log_L.append(float(split_line[1]))
                else:
                    disregard = True
            if 'e\t\tlog_L\t\tmaximised_overlap' in line:
                read_from_here = True
        if not disregard:
            cdf = rwt.cumulative_density_function(log_L)
            eccentricity_sample = rwt.pick_weighted_random_eccentricity(cdf, eccentricities)
            eccentricity_distribution.append(eccentricity_sample)
        else:
            print('disregarding file ' + single_file)

storage_file = args.results_folder + '/eccentricity_histogram_data.txt'
with open(storage_file, 'w') as f:
    for eccentricity in eccentricity_distribution:
        f.write(str(eccentricity) + '\n')

fig = plt.figure()
plt.hist(eccentricity_distribution, bins=np.logspace(-4, np.log10(0.2)))
plt.xlabel('eccentricity, $e$')

output_file = args.results_folder + '/eccentricity_histogram.pdf'
plt.savefig(output_file, bbox_inches='tight')
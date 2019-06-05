import python_code.split as split
import argparse


# Set up the argument parser
parser = argparse.ArgumentParser(
    "Split a result file into multiple smaller sets"
)
parser.add_argument('-r', '--result', help='Path of the result file to use')
parser.add_argument('-n', '--number-per-job', help='number of samples to analyse per job')
args = parser.parse_args()

# Split result fle into many subsets
result_file = args.result
number_per_job = int(args.number_per_job)
split.split_results_into_subsets(number_per_file=number_per_job, result_file=result_file)
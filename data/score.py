import argparse
import csv
import gzip
import math
import sys


def print_usage():
    usage = """NAME
       calc_score.py - Script to compare endpoint and global slack result.

SYNOPSIS
       calc_score.py OPTION... 

DESCRIPTION
       Calculate endpoint and global slack score with given reference result.

       Mandatory arguments:

       -r, --reference fname
              Name of the reference slack result file

       -t, --result fname
              Name of the output slack result file

       -s, --score score
              Full score for the current step

       Optional arguments:

       -e, --endpoints fname
              Name of the endpoint list file. If this option is provided, only vertices listed in this file are compared, otherwise all vertices in the reference file are compared.

       --threshold thres
              Value of the invalid result threshold

EXAMPLE
       python calc_score.py -r ref_slack.csv.gz -e end_list.list -t output_slack.csv.gz
"""
    print(usage)


def open_file(filename, mode):
    if mode == "read" or mode == "r":
        if filename.endswith('.gz'):
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')
    elif mode == "write" or mode == "w":
        if filename.endswith('.gz'):
            return gzip.open(filename, 'wt')
        else:
            return open(filename, 'w')
    else:
        raise ValueError(f"Unsupported mode {mode} in openFile")


def read_csv_file(fname):
    with open_file(fname, "r") as fh:
        reader = csv.reader(fh)
        data = {}
        skip_header = True
        for row in reader:
            if skip_header:
                skip_header = False
                continue
            data[row[0]] = [float(row[1]), float(row[2])]
        return data


def read_list_file(fname):
    with open_file(fname, "r") as fh:
        return [line.strip() for line in fh]


def calc_diff(result, ref):
    return (result - ref) / ref if ref != 0 else 0


def compare(ref_data, result_data, list_data, thres, full_score):
    compare_all = False
    if not list_data or len(list_data) == 0:
        compare_all = True

    max_abs_diff = -1e99
    max_abs_diff_pin = None
    all_diff = []
    total_cal_score = 0
    invalid_count = 0
    total_count = 0

    for pin, ref_slack in ref_data.items():
        total_count += 2
        result_slack = result_data.get(pin)

        if result_slack is None:
            invalid_count += 2
            continue

        for i in range(2):
            if math.isnan(result_slack[i]):
                invalid_count += 1
                continue

            diff = calc_diff(result_slack[i], ref_slack[i])
            abs_diff = abs(diff)

            if abs_diff > thres:
                invalid_count += 1
            else:
                total_cal_score += 1 - abs_diff

            if abs_diff > max_abs_diff:
                max_abs_diff = abs_diff
                max_abs_diff_pin = pin

            all_diff.append(diff)

    final_score = total_cal_score / total_count * full_score if total_count != 0 else 0
    invalid_pct = 100.0 * invalid_count / total_count if total_count != 0 else 0

    print(
        f"Total {total_count} slacks compared, step score: {final_score}, invalid slack count: {invalid_count}, invalid percent: {invalid_pct}")

    total_diff = sum(all_diff)
    diff_size = len(all_diff)
    diff_mean = total_diff / diff_size if diff_size != 0 else 0

    dev_total = sum([(diff - diff_mean) ** 2 for diff in all_diff])
    diff_stddev = math.sqrt(dev_total / diff_size) if diff_size != 0 else 0

    print(f"Max absolute slack difference pin: {max_abs_diff_pin}, abs slack difference: {max_abs_diff}")
    print(f"Slack difference distribution: mean = {diff_mean}, stddev = {diff_stddev}")


def main():
    parser = argparse.ArgumentParser(description='Script to compare endpoint and global slack result.')
    parser.add_argument('-r', '--reference', help='Name of the reference slack result file', required=True)
    parser.add_argument('-t', '--result', help='Name of the output slack result file', required=True)
    parser.add_argument('-s', '--score', help='Full score for the current step', type=int, required=True)
    parser.add_argument('-e', '--endpoints', help='Name of the endpoint list file')
    parser.add_argument('--threshold', help='Value of the invalid result threshold', type=float, default=0.2)
    args = parser.parse_args()

    ref_data = read_csv_file(args.reference)
    result_data = read_csv_file(args.result)
    list_data = read_list_file(args.endpoints) if args.endpoints else None

    compare(ref_data, result_data, list_data, args.threshold, args.score)


if __name__ == "__main__":
    main()

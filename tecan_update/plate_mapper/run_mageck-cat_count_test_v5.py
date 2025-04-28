import subprocess
import os
import pandas as pd
import glob
import csv
import argparse
import pathlib


# Python3
def mkargs(input_df):
    """
    Takes a list of all sample names and creates string format of labels
    and fastq files for mageck module. Labels will be separated by a
    ",", i.e. Day0_1,Day0_2. Fastq file string will have samples
    separated by a space, i.e. Day0_1 Day0_2 Day0_3
    """
    all_samples = input_df["Sample"].tolist()
    labels = str()
    fastq = ""
    count = 0
    for i in all_samples:
        i = str(
            i
        )  # Switches type from pandas unicode to python string. May not be necessary.
        if count == 0:
            labels += "{}".format(i)
            fastq += "fastq/" + i + ".fastq"
        else:
            labels += ",{}".format(i)
            fastq += " fastq/" + i + ".fastq"
        count += 1
    # print (labels)
    return labels, fastq


def cat_fastq_files(input_df):  # , GSF_dir
    # GSF_base_dir = "/research/dept/hart/PI_data_distribution/millergrp/GSF/"
    # src_dir = GSF_base_dir + GSF_dir
    # print(("copy from: {}".format(src_dir)))
    # copy_string = "find  " + src_dir + ' -name "*R1*.fastq.gz" -exec cp {} . \;'
    # print(copy_string)
    # os.system(str(copy_string))
    # os.system("gzip -d *.gz")
    print("\ninput_df:\n")
    print(input_df.head())
    print("\n  ***\n\n")
    fastq_dict = dict(list(zip(input_df.Sample, input_df.Index)))
    print(("fastq dict is {}".format(fastq_dict)))
    for key, value in fastq_dict.items():
        print(("key is {}, value is {}".format(key, value)))
        filecount = subprocess.check_output(
            f"ls -1 *{value}*.fastq | wc -l",
            shell=True,
            text=True,
            stderr=subprocess.STDOUT,
        ).strip()
        print(f"{filecount} fastq files found for {key}")
        if filecount == "1":
            # print("filecount == 1")
            os.system(f"mv *{value}*.fastq {key}.fastq")
        else:
            os.system(str("cat *{}*.fastq > {}.fastq".format(value, key)))

    return


def run_mageck_count(labels, fastq, lib_file):
    """
    Runs mageck count program from command line using os.system. Uses
    a string of labels, a string of fastq files that match to the labels,
    and the input file for mageck.
    Note: subprocess module does not work for more than one fastq file with spaces between.
    """
    # if dirc != None:
    if not os.path.exists("count"):
        os.mkdir("count")

    print("Running Mageck Count...\n")

    mageck_call = f'bsub -o count.log -P CAGE -R "rusage[mem=10000]" -q priority mageck count -l {lib_file} --sample-label {labels} --fastq {fastq} -n count/count'
    print(f" ***\n\n  MAGECK CALL:\n\n {mageck_call}\n\n ***\n\n")
    os.system(mageck_call)
    return


def run_mageck_test(input_df, read_count_file, dirc=None):
    """Runs mageck test program from linux command line. Only works for
    samples in triplicate with same prefix, i.e. "Day0" or "Day1"""

    # Changes directory to dirc argument. If user changes default, this
    #   function changes the current working directory

    curr_wd = os.getcwd()

    read_path = read_count_file
    read_filename = read_count_file.split("/")[-1]

    # print(f"{read_path}")
    # print(f"{read_filename}")

    if dirc != None:
        dirc = os.path.join(os.getcwd(), dirc)
        if not os.path.exists(dirc):
            os.mkdir(dirc)
        copy_command = f"cp {read_path} {dirc}/"
        os.system(copy_command)
        os.chdir(dirc)

    grouped_df = input_df.astype(str).groupby(["Set"])
    groups_list = []
    group_dict = {}
    for name, group in grouped_df:
        if str(name) == "Control":
            print("control")
            print(("print sample: {}".format(group["Sample"])))
            control_list = group["Sample"].tolist()
        else:
            print(("Name: {},  Group: {}".format(name, group)))
            sample_set_list = group["Sample"].tolist()
            groups_list.append(sample_set_list)
            group_dict[name] = sample_set_list

    for key, value in list(group_dict.items()):
        print(("Key: {}, value: {}".format(key, value)))
        mageck_test_call = f'bsub -o output.log -e error.log -P CAGE mageck test -k {read_filename} -t {",".join(group_dict[key])} -c {",".join(control_list)} -n {key} --pdf-report'

        os.system(mageck_test_call)
    os.chdir(curr_wd)
    return


def grab_input_files():
    pass


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process SRM library file")
    parser.add_argument("-i", help="Input key CSV file containing groups/samples")
    parser.add_argument(
        "-l",
        help="Input gRNA library file if running MaGeCK test",
        default="library.csv",
    )
    parser.add_argument(
        "-c",
        help='Choice "cat" to concatenate, "count" to count, "test" to test',
        type=str,
    )
    parser.add_argument("-k", help="gRNA read_count file", type=str)
    # parser.add_argument(
    #     "-d",
    #     help="destination directory (will create a new folder for tail of directory if one has not already been made)",
    #     type=str,
    #     default=os.getcwd(),
    # )
    parser.add_argument(
        "--GSF",
        help="Hartwell GSF directory to copy files from, just the part after /GSF",
        type=str,
        default="",
    )
    parser.print_help()
    args = parser.parse_args()
    lib_file = str(args.l)
    read_count_file = args.k
    in_file = args.i
    name_from_key = pathlib.PurePath(args.i).stem.split("key_")[1]
    input_df = pd.read_csv(in_file)
    print(("lib file is : {}".format(lib_file)))
    print(("directory is {}".format(name_from_key)))
    if args.c == "cat":
        cat_fastq_files(input_df)  # , args.GSF
    elif args.c == "count":
        labels, fastq = mkargs(input_df)
        samples = [i for i in input_df["Sample"].tolist() if pd.isnull(i) == False]
        run_mageck_count(labels, fastq, lib_file)
    elif args.c == "test":
        run_mageck_test(input_df, read_count_file, dirc=name_from_key)


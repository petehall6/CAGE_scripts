import os
import Bio
from Bio import SeqIO
import pandas as pd
import re, pdb
from io import StringIO
import numpy as np


def count_matches(df, summary_df, key):
    grouped = df.groupby("Name")
    match_dict = {}
    for grp in grouped:
        match_dict[grp[0]] = len(grp[1])

    summary_df[key] = pd.Series(match_dict)
    return summary_df


def _shift_index(x):
    # Helper function to format_sequences
    mismatch_list = x.split(",")
    new_list = []
    lappend = new_list.append
    for mismatch in mismatch_list:
        index, base = mismatch.split(":")
        new_index = int(index) + 1
        lappend(str(new_index) + ":" + base)

    return ",".join(new_list)


def format_sequences(temp_frame, output_fp):
    complement_dict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "V": "V"}

    temp_frame["mismatches"] = temp_frame["mismatches"].apply(_shift_index)

    rev_strand_index = temp_frame["strand"] == "-"

    temp_frame["sequence"].update(
        temp_frame["sequence"]
        .loc[rev_strand_index]
        .apply(lambda x: "".join([complement_dict[b] for b in x[::-1]]))
    )

    temp_frame.to_csv(output_fp, sep="\t", index=False)
    return None


def count_alignment_and_format_bowtie(summary_df, bowtie_fp_list, no_format=False):
    """
    Alters the passed summary_df inplace with count info of passed bowtie aligments
    Also formats the bowtie output and writes to another file
    """
    BOWTIE_FORMAT_APPEND_STR = ".formatted"

    for bowtie_fp in bowtie_fp_list:
        match_type = re.search(r"CRISPR_([A-Za-z]+_\d)", bowtie_fp).group(1)

        alignment_df = pd.read_csv(
            bowtie_fp,
            delimiter="\t",
            names=["Name", "strand", "chrom", "pos", "sequence", "mismatches"],
            header=None,
        )

        summary_df = count_matches(alignment_df, summary_df, match_type)

        # shift the mismatch index to start from 1 and reverse complement
        # - strand matches
        if not no_format:
            format_sequences(alignment_df, bowtie_fp + BOWTIE_FORMAT_APPEND_STR)

    return summary_df


def call_bowtie_for_OTA(OTA_choice, for_primers=False):
    # look for mismatches from the full 23 bp sequence
    # summarize the number of hits of each type.
    # make summary frame
    bowtie_species_paths = {
        "human": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/GCA_000001405.15_GRCh38_no_alt_analysis_set",
        "mouse": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/mm10",
        "rat": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/rn5",
        "zebrafish": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/Zv9",
        "hamster": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/CHO-K1_2012",
        "horse": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/EquCab2.0",
        "monkey": "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/chlorocebus_sabaeus/Csabaeus_Vero",
    }

    bowtie_args_sets = [
        {
            "v_arg": 1,
            "nsub": "CRISPR_nsub.faa",
            "final_bowtie": "CRISPR_long_0mm.bowtie",
            "easy_name": "long_0",
        },
        {
            "v_arg": 2,
            "nsub": "CRISPR_nsub.faa",
            "final_bowtie": "CRISPR_long_1mm.bowtie",
            "easy_name": "long_1",
        },
        {
            "v_arg": 3,
            "nsub": "CRISPR_nsub.faa",
            "final_bowtie": "CRISPR_long_2mm.bowtie",
            "easy_name": "long_2",
        },
        {
            "v_arg": 1,
            "nsub": "CRISPR_nsub_14.faa",
            "final_bowtie": "CRISPR_short_0mm.bowtie",
            "easy_name": "short_0",
        },
    ]

    if for_primers:
        bowtie_args_sets = [
            {
                "v_arg": 1,
                "nsub": "primers.fa >",
                "final_bowtie": "primer_0mm.bowtie",
                "easy_name": "0_mm",
            },
            {
                "v_arg": 2,
                "nsub": "primers.fa >",
                "final_bowtie": "primer_1mm.bowtie",
                "easy_name": "1_mm",
            },
            {
                "v_arg": 3,
                "nsub": "primers.fa >",
                "final_bowtie": "primer_2mm.bowtie",
                "easy_name": "2_mm",
            },
        ]

    if OTA_choice in bowtie_species_paths:
        bowtie_path = bowtie_species_paths[OTA_choice]
        for set in bowtie_args_sets:
            bowtie_call = f"bowtie  --offbase 1 -v {set['v_arg']}  -a -y --best --suppress 6,7 -f {bowtie_path} {set['nsub']} {set['final_bowtie']}"
            # print(f"\nBowtie Call:\n   * {bowtie_call}\n")
            print(f"\n{set['easy_name']} bowtie call:")
            os.system(bowtie_call)
        print()
        return None
    else:
        print("Not ready yet.")
        return None


def off_target_all(OTA_choice, no_format=False, casType="cas9"):
    CAS_TYPES = ["cas9", "cas12a"]
    if casType not in CAS_TYPES:
        print(f"Incorrect cas type supplied to off_target_all: {casType}")
        print("Please correct and re-run the program.")
        return 1

    inputfilename = "CRISPR.fa"

    complement_dict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    crispr_sequences = list(SeqIO.parse(inputfilename, "fasta"))

    # make summary frame

    summary_index = list()
    columns = ["gRNA", "long_0", "long_1", "long_2", "short_0"]
    for seq_record in crispr_sequences:
        summary_index.append(seq_record.name)

    summary_frame = pd.DataFrame(columns=columns, index=summary_index)
    summary_frame = summary_frame.fillna(0)
    # add sequences to summary frame
    for seq_record in crispr_sequences:
        summary_frame.loc[seq_record.name, "gRNA"] = str(seq_record.seq)

    SeqIO.write(crispr_sequences, "CRISPR_nsub_14.faa", "fasta")
    for seq_record in crispr_sequences:
        seq_record.seq = seq_record.seq.tomutable()
        if casType == "cas9":
            seq_record.seq[20] = "N"
        if casType == "cas12a":
            seq_record.seq[3] = "N"

    SeqIO.write(crispr_sequences, "CRISPR_nsub.faa", "fasta")

    for seq_record in crispr_sequences:
        if casType == "cas9":
            seq_record.seq = seq_record.seq[8:]
        if casType == "cas12a":
            seq_record.seq = seq_record.seq[0:16]

    SeqIO.write(crispr_sequences, "CRISPR_nsub_14.faa", "fasta")

    # look for mismatches from the full 23 bp sequence
    # summarize the number of hits of each type.
    # make summary frame
    call_bowtie_for_OTA(OTA_choice)

    summary_frame = count_alignment_and_format_bowtie(
        summary_frame,
        [
            "CRISPR_short_0mm.bowtie",
            "CRISPR_long_0mm.bowtie",
            "CRISPR_long_1mm.bowtie",
            "CRISPR_long_2mm.bowtie",
        ],
        no_format=no_format,
    )

    # Convert OTA results to int from float
    col_dtypes_dict = {
        "long_0": "Int64",
        "long_1": "Int64",
        "long_2": "Int64",
        "short_0": "Int64",
    }
    summary_frame = summary_frame.astype(col_dtypes_dict)

    summary_frame = summary_frame.sort_values(["long_0", "long_1", "long_2", "short_0"])
    return summary_frame


def off_target_all_Cas12a(OTA_choice, no_format=False):
    return off_target_all(OTA_choice, no_format, casType="cas12a")


def off_target_all_primers(directory, OTA_choice, no_format=False):
    # TODO: Investigate removal of "directory" argument.
    inputfilename = "primers.fa"
    # directory='/rgs01/project_space/millergrp/CrisprDesigns/common/local/crispr_design/python3/'
    complement_dict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    with open(inputfilename) as f:
        lines = f.read().splitlines()
    lines = [item.strip(">") for item in lines]
    print("\nPrimer off targets:")
    column_list = ["primer_name", "strand", "chr", "start", "sequence"]
    summary_cols = ["primer_name", "sequence"]
    primer_dict = dict(zip(lines[::2], lines[1::2]))
    summary_frame = pd.DataFrame(list(primer_dict.items()), columns=summary_cols)
    call_bowtie_for_OTA(OTA_choice, for_primers=True)

    primer_OTA_list = ["primer_0mm.bowtie", "primer_1mm.bowtie", "primer_2mm.bowtie"]
    for index, OTA in enumerate(primer_OTA_list):
        the_file = os.path.join(directory, OTA)
        temp_df = pd.read_csv(the_file, index_col=False, sep="\t", names=column_list)
        temp_df = (
            temp_df.groupby(["primer_name"])
            .size()
            .reset_index(name=str(str(index) + "_mm"))
        )
        summary_frame = pd.merge(summary_frame, temp_df, on="primer_name")
    print(summary_frame)
    summary_frame.to_csv("primer_OTA.txt", sep="\t", index=False)
    return summary_frame
#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Literal

import pandas as pd
import janitor
import pyranges as pr

__version__ = "0.0.1"

# TODO: Standardize column names
# TODO: In case of a segment spanning the arms, what to do? It can cause
# problems if the sum cn become more than 100% of the genome, perhaps
# include the normalized length once or just do over the whole chromosome

def get_chromosomal_arm_lengths(genome: str="hg38") -> pd.DataFrame:

    url = f"http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/cytoBand.txt.gz"
    table = (
        pd.read_table(url, names=["chrom","chromStart","chromEnd",
                                  "arm","gieStain"] )
        .transform_column("arm", lambda x: x.str[0], elementwise=False)
        .dropna(subset=["arm"])
        .groupby(["chrom", "arm"], as_index=False )
        .agg({"chromStart": "min", "chromEnd": "max"})
        .sort_naturally("chrom")
        .rename(columns={"chrom": "Chromosome", "chromStart": "Start",
                         "chromEnd": "End", "arm": "Name"})
    )

    table = table.assign(ArmLength=table.apply(lambda x: x.End - x.Start, axis=1))

    return table


def ucsc_to_ncbi(dataframe: pd.DataFrame) -> pd.DataFrame:

    result = dataframe.transform_column("chromosome",
                                        lambda x: x.str.replace, "chr", "",
                                        elementwise=False)

    return result

def ncbi_to_ucsc(dataframe: pd.DataFrame) -> pd.DataFrame:

    result = (
        dataframe
        .astype({"chromosome": "str"})
        .transform_column("chromosome", lambda x: x.radd("chr"),
                          elementwise=False)
    )

    return result


def classifier(row) -> Literal['S', 'U', 'HU'] | None:
    cnb = row.cnb

    if pd.isnull(cnb): # No CNB calculated, we assume stable
        return "S"

    length = row.normalized_length
    if length > 0.95 and cnb < 2.5:
        return "S"
    if length > 0.95 and cnb > 2.5:
        return "U"
    if length < 0.95 and cnb > 27:
        return "HU"
    if length < 0.95 and cnb < 27:
        return "U"


def copy_number_burden(value) -> float:
    cnb = value.sum() / 3e9
    cnb = round(cnb * 100, 4)
    return cnb


def adjust_call(value: int, ploidy: int) -> float:

    cn = value["absolute_cn"]
    adjusted = cn - ploidy
    return adjusted


def harmonize_columns(
    dataframe: pd.DataFrame,
        format: Literal['ascat_sc', 'ichorcna', 'ace']) -> pd.DataFrame:

    match format:
        case "ascat_sc":
            dataframe = dataframe.rename(
                columns={"total_copy_number": "absolute_copy_number"}
            )
        case "ace":
            dataframe = dataframe.rename(
                columns={"copies": "absolute_copy_number"}
            )
        case "ichorcna":
            dataframe = dataframe.rename(
                columns={"copy.number": "absolute_copy_number",
                         "chr": "chromosome"}
            )
        case _:
            raise ValueError(f"Unsupported file format {format}")

    return dataframe

def normalized_scna_length(absolute_calls: pd.DataFrame,
                           arm_data: pd.DataFrame) -> pd.DataFrame:
    # Convert the dataframes to PyRanges objects
    calls_pr = pr.PyRanges(absolute_calls.rename(
        columns={"chromosome": "Chromosome", "start": "Start", "end": "End"}))
    arms_pr = pr.PyRanges(arm_data.rename(
        columns={"chrom": "Chromosome", "chromStart": "Start",
                 "chromEnd": "End", "arm": "Name"}))

    # Perform intersection
    intersected = calls_pr.join(arms_pr, how='left')

    # Calculate total chromosome lengths from arm_data
    chromosome_lengths = arm_data.groupby('Chromosome', as_index=False).agg(
        TotalLength=('End', 'max'))
    chromosome_lengths.rename(columns={'chrom': 'Chromosome'}, inplace=True)

    # Calculate normalized length
    intersected_df = intersected.df
    intersected_df['new_start'] = intersected_df[['Start', 'Start_b']].max(
        axis=1)
    intersected_df['new_end'] = intersected_df[['End', 'End_b']].min(axis=1)
    intersected_df['length'] = (intersected_df['new_end'] -
                                intersected_df['new_start'])
    intersected_df['normalized_length'] = intersected_df['length'] / intersected_df['ArmLength']

    # Merge with chromosome lengths to get the total length for normalization
    intersected_df = pd.merge(intersected_df, chromosome_lengths,
                             on='Chromosome', how='left')

    # Normalization: check if the segment spans both arms. If so,
    # normalize over the total chromosome length.
    intersected_df['normalized_length'] = intersected_df.apply(
        lambda x: x['length'] / x['TotalLength']
       if (x['Start'] <= x['Start_b'] and x['End'] >= x['End_b'])
       else x['length'] / x['ArmLength'], axis=1)

    # Clean up the DataFrame using pyjanitor to standardize column names
    cleaned_df = (
        intersected_df
        .rename(columns={"Name": "arm", "ArmLength": "arm_length"})
        .select_columns(['Chromosome', 'Start', 'End', 'sample', 'arm',
                         'arm_length', 'normalized_length'])
        .clean_names()
    )

    return cleaned_df


def call_patterns(segments: pd.DataFrame, arm_data: pd.DataFrame,
                         ploidy: int):

    segment_length = (
        normalized_scna_length(segments, arm_data)
        .merge(segments, on=["chromosome", "start", "end", "sample"],
               how="inner")
    )

    segment_length = (
        segment_length
        .assign(call=segment_length.apply(adjust_call, ploidy=ploidy, axis=1))
    )


    normalized_length = (
        segment_length
        .groupby("sample")["normalized_length"]
        .quantile(0.75).to_frame()
    )

    cnb = (
        segment_length
        .query("call != 0")
        .groupby("sample")
        .agg({"length": copy_number_burden})
        .rename(columns={"length": "cnb"})
    )

    merged_table = pd.concat([normalized_length, cnb], axis=1)

    merged_table = merged_table.assign(
        pattern=merged_table.apply(classifier, axis=1)
    )

    return merged_table


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--ploidy", type=int, help="Ploidy of the sample",
                        default=2)
    parser.add_argument("--genome", choices=("hg19", "hg38"), default="hg38")
    parser.add_argument("--file-format", choices=("ascat_sc", "ichorcna",
                                                  "ace"))
    parser.add_argument("--genome-style", choices=("ucsc", "ncbi"),
                        default="ucsc")
    parser.add_argument("--sample-name", default="Sample")
    parser.add_argument("source", help="Source segment file")
    parser.add_argument("destination", help="Destination to save to")
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)

    options = parser.parse_args()

    source = pd.read_table(options.source).clean_names()

    if options.sample:
        source["sample"] = options.sample
    # Fall back to integrated sample, except for ACE that doesn't have it
    elif options.file_format == "ace":
        source["sample"] = Path(options.source).stem

    source = harmonize_columns(source, options.file_format).astype(
        {"chromosome": "str"}
    )

    if options.genome_style == "ucsc" and not source.chromosome.str.startswith("chr").all():
        source = ncbi_to_ucsc(source)

    if options.genome_style == "ncbi"and not source.chromosome.str.startswith("chr").any():
        source = ucsc_to_ncbi(source)

    source = source.assign(length=source.end.sub(source.start))

    destination = options.destination

    arm_data = get_chromosomal_arm_lengths(genome=options.genome)

    classified_sample = call_patterns(source, arm_data,
                                      options.ploidy)

    classified_sample.to_csv(f"{destination}.txt", sep="\t",
                             index=True)


if __name__ == "__main__":
    main()

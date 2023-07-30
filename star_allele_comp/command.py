"""
This module will be called from command line.
It extracts the data from csv, compare the data and print/save the results.
"""
import argparse
from collections import defaultdict
from typing import Iterable

import pandas as pd

from .comparison import compare_method
from .summary import (
    save_all_summary,
    print_all_summary,
)
from .plot import plot_summary


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser"""
    parser = argparse.ArgumentParser(description="Star Allele Comparator")
    parser.add_argument(
        "csv",
        nargs="+",
        metavar="file1.csv file2.csv ...",
        help="CSV files for allele calling data",
    )
    parser.add_argument(
        "--ref",
        metavar="ground_truth_method_name",
        help="Reference method for comparison (If not specified, reference method is the first method)",
    )
    parser.add_argument(
        "--family",
        choices=["HLA", "KIR", "hla", "kir"],
        default="HLA",
        help="Choose 'HLA' or 'KIR' family for comparison",
    )

    # Optional flags
    # parser.add_argument("--summary-type", nargs="+", help="Restrict summary type to print/save")
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument(
        "--save",
        nargs="?",
        metavar="path/to/save",
        help="Saving all results in [path].xxx.oo",
        const="output",  # the default if '--save [empty]'
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable output to screen"
    )
    return parser


def read_file(filename: str) -> Iterable[tuple[str, str, list[str]]]:
    """Read data from csv"""
    if filename.endswith(".tsv"):
        df_raw = pd.read_csv(filename, sep="\t")
    else:
        df_raw = pd.read_csv(filename)
    columns = [str(i).lower() for i in df_raw.columns]
    df_raw = df_raw.set_axis(columns, axis=1)
    column_allele = [i for i in columns if i.startswith("allele") and i != "alleles"]

    if "method" not in df_raw.columns:
        df_raw["method"] = filename

    df_raw = df_raw[["method", "id", *column_allele]]

    for _, row in df_raw.iterrows():
        yield row.method, row.id, [
            row[col].strip()
            for col in column_allele
            if not pd.isnull(row[col]) and row[col].strip()
        ]
        if "alleles" in df_raw.columns:
            yield row.method, row.id, [
                s.strip() for s in row.alleles.split("_") if s.strip()
            ]


def files_to_cohort(files: Iterable[str]) -> dict[str, dict[str, list[str]]]:
    """Turn multiple csv files into CohortInput."""
    cohort: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for filename in files:
        for method, sample_id, alleles in read_file(filename):
            cohort[method][sample_id] += [
                col.strip() for col in alleles if not pd.isnull(col)
            ]
    return cohort


def entrypoint(commands: list[str] | None = None) -> None:
    """The function will trigger by command line"""
    parser = create_parser()
    args = parser.parse_args(commands)

    cohort = files_to_cohort(args.csv)
    if args.ref:
        ref_method = args.ref
    else:
        ref_method = list(cohort.keys())[0]

    result = compare_method(cohort, ref_method, args.family)
    result_df = result.to_dataframe()

    path_save = args.save
    if path_save and path_save.endswith("/"):
        path_save += "output"

    if args.verbose:
        print(result)
        print_all_summary(result_df)

    if path_save:
        with open(path_save + ".summary.txt", "w") as f_handle:
            print_all_summary(result_df, file_handler=f_handle)
            print(result, file=f_handle)
        save_all_summary(result_df, path_save)

    if args.plot and path_save:
        figs = plot_summary(result_df)
        for i, fig in enumerate(figs):
            fig.write_image(path_save + f".figure.{i}.png")


if __name__ == "__main__":
    entrypoint()

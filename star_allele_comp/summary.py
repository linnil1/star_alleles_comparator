"""
Summary functions for the comparison result from comparison.py.
"""
from __future__ import annotations
from typing import Any, Callable, TextIO, Iterable
from functools import wraps
import sys

import numpy as np
import pandas as pd


def apply_method_wise(
    func: Callable[..., pd.DataFrame],
) -> Callable[..., pd.DataFrame]:
    """
    A decorator that applies the given function to each method
    if the 'method' column exists in the DataFrame.
    """

    @wraps(func)
    def wrapper(df_any: pd.DataFrame, **kwargs: Any) -> pd.DataFrame:
        if "method" in df_any.columns and len(set(df_any["method"])) > 1:
            return df_any.groupby(["method"]).apply(
                lambda df_sub: func(df_sub, **kwargs)
            )
        # not need to group by method
        return func(df_any, **kwargs)

    return wrapper


@apply_method_wise
def table_confusion(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get Confusion Matrix

    Returns:
    ```
    MatchRes    1  2  3 FN FP
    RefRes
    -1          0  0  0  0  1
    1           0  0  0  1  0
    2           2  6  0  1  0
    3           0  0  1  0  0
    4           0  0  1  0  0
    ```
    """
    confusion_df = results_df.copy()
    confusion_df["Count"] = 1
    confusion_df = confusion_df.pivot_table(
        index=["allele1_res"],
        columns=["match_res"],
        values=["Count"],  # type: ignore
        aggfunc=sum,
        fill_value=0,
    )
    confusion_df = confusion_df.astype(pd.Int64Dtype())
    confusion_df.index.name = "ref_res"
    confusion_df = confusion_df.droplevel(0, axis=1)
    return confusion_df


@apply_method_wise
def table_summarize(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get summary of the accuracy in each resolution

    * The 'Resolution' column represents the resolution level for each comparison.
    * 'num_match' represents the number of alleles that match the alleles in the ground truth method
        under the specific resolution level.
    * 'num_ref' indicates the number of reference alleles with a resolution level greater than or equal to the resolution level in the comparison.
    * 'Accuracy' is calculated as the ratio of 'num_match' to 'num_ref'.
    * 'Accuracy' in FP represents the False Discovery Rate (FDR).
    * 'Accuracy' in FN represents the False Negative Rate (FNR).

    Returns:
    ```
        Resolution  num_match  num_ref  Accuracy
    0          4          0        2  0.000000
    1          3          0        2  0.000000
    2          2          0        2  0.000000
    3          1          2        2  1.000000
    4          0          2        2  1.000000
    5         FP          2        0  0.333333
    6         FN          0        0  0.000000
    ```
    """
    max_res = max(*results_df["allele1_res"], *results_df["match_res"])
    summary_list: list[tuple[Any, ...]] = []
    for resolution in range(max_res, -1, -1):
        res_ref = len(results_df[results_df["allele1_res"] >= resolution])
        res_match = len(results_df[results_df["match_res"] >= resolution])
        summary_list.append(
            (resolution, res_match, res_ref, np.divide(res_match, res_ref))
        )

    n_fp = len(results_df[results_df["match_str"] == "FP"])
    n_fn = len(results_df[results_df["match_str"] == "FN"])
    # FP -> FDR False Discover Rate
    # FN -> FNR False Negative Rate
    summary_list.extend(
        [
            (
                "FP",
                n_fp,
                0,
                np.divide(n_fp, n_fp + len(results_df)),
            ),
            (
                "FN",
                n_fn,
                0,
                np.divide(n_fn, n_fn + len(results_df)),
            ),
        ]
    )

    df_res = pd.DataFrame(
        summary_list, columns=["Resolution", "num_match", "num_ref", "Accuracy"]
    )
    df_res["num_match"] = df_res["num_match"].fillna(0).astype(pd.Int64Dtype())
    df_res["num_ref"] = df_res["num_ref"].fillna(0).astype(pd.Int64Dtype())
    df_res["Resolution"] = df_res["Resolution"].astype(str)
    # loop df_res and print
    # for index, row in df_res.iterrows():
    #     if row["Resolution"] == "FP":
    #         print(f"FP: {row['num_match']:3d}")
    #         continue
    #     print(
    #         f"{row['Resolution']:2d}: {row['num_match']:3d} / {row['num_ref']:3d} = {row['Accuracy']:.3f}"
    #     )
    # print(f"{resolution:2d}: {res_match:3d} / {res_ref:3d} = {np.divide(res_match,res_ref):.3f}")
    # print(f"FP: {:3d}")
    return df_res


@apply_method_wise
def table_summarize_group_by(
    results_df: pd.DataFrame, group_by: str | list[str]
) -> pd.DataFrame:
    """
    This function is similar to table_summarize
    but operates on a grouped DataFrame by `group_by` column.
    """
    results_summary_df = results_df.groupby(group_by).apply(table_summarize)
    results_summary_df = results_summary_df.droplevel(1).reset_index()
    return results_summary_df


def compact_summary(df_res: pd.DataFrame, group_by: Iterable[str] = []) -> pd.DataFrame:
    """
    Transform summary in table_summarize into one line.

    Examples:
    ```
        Resolution  num_match  num_ref  Accuracy
    0          3          2        2  1.000000
    1          2          8       11  0.727273
    2          1         10       13  0.769231
    3          0         11       13  0.846154
    4         FP          2        0       NaN
    5         FN          2        0       NaN

    Transform into:

    Accuracy                          num_match                num_ref
    0         1         2         3    0   1   2  3  FN FP     0   1   2  3 FN FP
    0.846154  0.769231  0.727273  1.0  11  10  8  2  2  2      13  13  11 2 0  0
    ```
    """
    if "method" in df_res.columns or "method" in df_res.index.names:
        remove_index = False
    else:
        remove_index = True
        df_res["method"] = "dummy"

    df_res_one_line: pd.DataFrame = df_res.pivot_table(
        index=["method", *group_by],
        columns=["Resolution"],
        values=["num_match", "Accuracy", "num_ref"],  # type: ignore
        aggfunc="first",
    )
    df_res_one_line["num_match"] = df_res_one_line["num_match"].fillna(0)
    df_res_one_line["num_ref"] = df_res_one_line["num_ref"].fillna(0)
    if remove_index:
        if df_res_one_line.index.nlevels > 1:
            df_res_one_line = df_res_one_line.droplevel(0)
        else:
            # df_res_one_line = df_res_one_line.iloc[0]
            df_res_one_line = df_res_one_line.reset_index(drop=True)
    return df_res_one_line


def get_summary_meta(summary_type: Iterable[str] = []) -> dict[str, dict[str, Any]]:
    """Get title, corresponding functions for each summary method"""

    def idenity(i: Any) -> Any:
        return i

    meta = {
        "raw": {
            "title": "# Raw dataframe",
            "func": idenity,
            "compact": idenity,
            # = df
        },
        "accuracy": {
            "title": "# Resolution Accuracy summary (one-line version)",
            "func": table_summarize,
            "compact": compact_summary,
            # = table_summarize(df)
            # = compact_summary(table_summarize(df))
        },
        "confusion": {
            "title": "# Confusion matrix",
            "func": table_confusion,
            "compact": idenity,
            # = table_confusion(df)
        },
        "accuracy-sample": {
            "title": "# Sample Accuracy summary per resolution",
            "func": lambda df: table_summarize_group_by(df, group_by="id"),
            "compact": lambda df: compact_summary(df, group_by=["id"]),
            # = table_summarize_group_by(df, group_by="id")
            # = compact_summary(table_summarize_group_by(df, group_by="id"), group_by=["id"])
        },
        "accuracy-gene": {
            "title": "# Gene Accuracy summary per resolution",
            "func": lambda df: table_summarize_group_by(df, group_by="gene"),
            "compact": lambda df: compact_summary(df, group_by=["gene"]),
            # = table_summarize_group_by(df, group_by="gene")
            # = compact_summary(table_summarize_group_by(df, group_by="gene"), group_by=["gene"])
        },
    }
    if not summary_type:
        return meta
    return {st: meta[st] for st in summary_type if st in meta}


def print_all_summary(
    df_result: pd.DataFrame,
    include_raw: bool = False,
    file_handler: TextIO = sys.stdout,
) -> None:
    """
    Print all the summary in the dataframe.

    Args:
        df_result: The DataFrame contained from `compare_cohort` or `compare_method`.
        include_raw: If True, print the raw dataframe of allele-by-allele comparison.
        file_handler: Output to stdout by default
    """
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", 60)
    pd.set_option("display.precision", 3)

    print("# Samples: ", len(df_result["id"].unique()), file=file_handler)
    print(df_result["id"].unique(), file=file_handler)
    if "method" in df_result.columns:
        print("# Methods:", len(df_result["method"].unique()), file=file_handler)
        print(df_result["method"].unique(), file=file_handler)
    print(file=file_handler)

    for metric_name, meta in get_summary_meta().items():
        if not include_raw and "raw" == metric_name:
            continue
        print(meta["title"], file=file_handler)
        print(meta["compact"](meta["func"](df_result)), file=file_handler)
        print(file=file_handler)


def save_all_summary(df_results: pd.DataFrame, path: str, compact: bool = True) -> None:
    """
    Save all summary to CSV files.

    Args:
        df_results: The DataFrame containing the summary data.
        path: The prefix of the file path. The summary will be saved to path.xx.oo.csv.
        compact: If True, save the compact version of the summary with the suffix ".compact.csv".
                 Defaults to True.
    """
    for metric_name, meta in get_summary_meta().items():
        df_summary = meta["func"](df_results)
        df_summary.to_csv(f"{path}.{metric_name}.csv")
        if not compact:
            continue
        df_compact = meta["compact"](df_summary)
        if len(df_compact.columns.names) > 1:
            df_compact.columns = ["-".join(map(str, col)) for col in df_compact.columns]
            df_compact.to_csv(f"{path}.{metric_name}.compact.csv")

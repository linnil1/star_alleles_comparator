"""
This module contains functions to plot results.
"""
from typing import Iterable
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .summary import table_summarize, table_summarize_group_by


def plot_bar(summary_df: pd.DataFrame, x: str, color: str | None) -> go.Figure:
    """Create a bar plot of summary dataframe"""
    summary_df = summary_df.reset_index()
    # print(summary_df#)
    category_orders = {
        "Resolution": sorted(summary_df["Resolution"].unique()),
    }
    if "method" in summary_df.columns:
        category_orders["method"] = sorted(summary_df["method"].unique())
    if "gene" in summary_df.columns:
        category_orders["gene"] = sorted(summary_df["gene"].unique())
    if "id" in summary_df.columns:
        category_orders["id"] = sorted(summary_df["id"].unique())

    return (
        px.bar(
            summary_df,
            x=x,
            y="Accuracy",
            color=color,
            text="Accuracy",
            text_auto=".2f",
            barmode="group",
            category_orders=category_orders,
        )
        .update_layout(
            yaxis_title="Accuracy/TNR/TPR",
        )
        .update_traces(
            textposition="outside", marker_line_width=1, marker_line_color="black"
        )
    )


def plot_summary(results_df: pd.DataFrame) -> Iterable[go.Figure]:
    """
    Generate multiple bar plots based on the result DataFrame obtained from compare_cohort or compare_method.
    """
    summary_df = table_summarize(results_df)
    if not "method" in results_df.columns:
        yield plot_bar(summary_df, x="Resolution", color=None).update_layout(
            title="Resolution Accuracy"
        )
    else:
        yield plot_bar(summary_df, x="Resolution", color="method").update_layout(
            title="Resoution Accuracy"
        )

    summary_gene_df = table_summarize_group_by(results_df, group_by=["gene"])
    if not "method" in results_df.columns:
        yield plot_bar(summary_gene_df, x="gene", color="Resolution").update_layout(
            title="Gene Accuracy"
        )
        # yield plot_bar(summary_gene_df, x="Resolution", color="gene")
    else:
        resolution_value = sorted(summary_gene_df["Resolution"].unique())
        for res in resolution_value:
            yield plot_bar(
                summary_gene_df[summary_gene_df["Resolution"] == res],
                x="gene",
                color="method",
            ).update_layout(title=f"Gene Accuracy (Resolution = {res})")

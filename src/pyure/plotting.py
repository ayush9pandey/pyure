"""Plotting helpers for simulation output."""

from __future__ import annotations

import re

import pandas as pd
import plotly.graph_objects as go


def parse_species_input(raw_text: str) -> list[str]:
    """Parse comma-separated or line-separated species names."""
    if not raw_text:
        return []
    names = [
        item.strip()
        for item in re.split(r"[,\n\r]+", raw_text)
        if item.strip()
    ]
    return list(dict.fromkeys(names))


def validate_requested_species(
    requested: list[str],
    available_columns: list[str],
) -> tuple[list[str], list[str]]:
    """Return found and missing species names."""
    available = set(available_columns)
    found = [name for name in requested if name in available]
    missing = [name for name in requested if name not in available]
    return found, missing


def default_species_column(dataframe: pd.DataFrame) -> str | None:
    """Return the first non-time dataframe column."""
    for column in dataframe.columns:
        if str(column).lower() != "time":
            return str(column)
    return None


def build_species_plot(dataframe: pd.DataFrame, species: list[str]) -> go.Figure:
    """Build a clean multi-trace time series figure."""
    x_column = "time" if "time" in dataframe.columns else dataframe.columns[0]
    figure = go.Figure()
    for name in species:
        figure.add_trace(
            go.Scatter(
                x=dataframe[x_column],
                y=dataframe[name],
                mode="lines",
                name=name,
                line={"width": 2},
            )
        )

    figure.update_layout(
        template="plotly_white",
        height=460,
        margin={"l": 48, "r": 24, "t": 36, "b": 48},
        xaxis_title="Time",
        yaxis_title="Concentration",
        legend_title_text="Species",
        hovermode="x unified",
    )
    return figure

"""Initial-condition helpers for generated models and simulations."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

DEFAULT_INITIAL_CONDITION_CSV = "PURE_TXTL_initial_values_Final.csv"


def load_initial_conditions(path: Path | None = None) -> dict[str, float]:
    """Load non-zero default initial conditions from the PURE notebook CSV."""
    csv_path = path or Path(__file__).resolve().parents[2] / DEFAULT_INITIAL_CONDITION_CSV
    if not csv_path.exists():
        return {
            "T7RNAP": 1.0,
            "DNA": 0.005,
            "ATP": 3750.0,
            "GTP": 2500.0,
            "CTP": 1250.0,
            "UTP": 1250.0,
        }

    initial_conditions: dict[str, float] = {}
    with csv_path.open(mode="r", newline="", encoding="utf-8-sig") as infile:
        reader = csv.reader(infile)
        for row in reader:
            if len(row) < 2:
                continue
            species = row[0].strip()
            try:
                value = float(row[1])
            except ValueError:
                continue
            if species and value != 0:
                initial_conditions[species] = value
    return initial_conditions


def default_initial_conditions_for_species(
    species_names: list[str],
    fallback_dna_species: str | None = None,
) -> dict[str, float]:
    """Return defaults limited to species that exist in the current model."""
    available = set(species_names)
    defaults = {
        species: value
        for species, value in load_initial_conditions().items()
        if species in available
    }

    if fallback_dna_species and fallback_dna_species in available:
        defaults.setdefault(fallback_dna_species, 0.005)

    return dict(sorted(defaults.items()))


def initial_conditions_to_dataframe(initial_conditions: dict[str, float]) -> pd.DataFrame:
    """Convert an initial-condition dict to an editable dataframe."""
    return pd.DataFrame(
        [
            {"Species": species, "Initial concentration": float(value)}
            for species, value in sorted(initial_conditions.items())
        ]
    )


def dataframe_to_initial_conditions(dataframe: pd.DataFrame) -> dict[str, float]:
    """Convert an edited dataframe back into a simulation initial-condition dict."""
    initial_conditions: dict[str, float] = {}
    if dataframe is None or dataframe.empty:
        return initial_conditions

    for _, row in dataframe.iterrows():
        species = str(row.get("Species", "")).strip()
        if not species:
            continue
        value = float(row.get("Initial concentration", 0.0))
        if value != 0:
            initial_conditions[species] = value
    return initial_conditions

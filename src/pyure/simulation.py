"""Bioscrape simulation wrapper."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pyure.exceptions import SimulationError


@dataclass(frozen=True)
class SimulationSettings:
    start: float = 0.0
    end: float = 3600.0
    points: int = 200
    stochastic: bool = False
    safe: bool = False


def simulate_sbml(
    sbml_xml: str,
    settings: SimulationSettings | None = None,
    initial_conditions: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Run a Bioscrape simulation from SBML XML and return a dataframe."""
    settings = settings or SimulationSettings()
    if settings.end <= settings.start:
        raise SimulationError("Simulation end time must be greater than start time.")
    if settings.points < 2:
        raise SimulationError("Simulation needs at least two time points.")

    try:
        from bioscrape.sbmlutil import read_model_from_sbml
        from bioscrape.simulator import py_simulate_model
    except Exception as exc:
        raise SimulationError(
            "Bioscrape is not available in this Python environment."
        ) from exc

    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False, mode="w", encoding="utf-8") as handle:
        handle.write(sbml_xml)
        sbml_path = Path(handle.name)

    try:
        model = read_model_from_sbml(str(sbml_path))
        if initial_conditions:
            available_species = set(model.get_species_list())
            filtered_initial_conditions = {
                species: float(value)
                for species, value in initial_conditions.items()
                if species in available_species
            }
            if filtered_initial_conditions:
                model.set_species(filtered_initial_conditions)
        timepoints = np.linspace(settings.start, settings.end, settings.points)
        result = py_simulate_model(
            timepoints,
            Model=model,
            stochastic=settings.stochastic,
            safe=settings.safe,
            return_dataframe=True,
        )
        if not isinstance(result, pd.DataFrame):
            result = pd.DataFrame(result)
        if "time" not in result.columns:
            result.insert(0, "time", timepoints)
        elif result.columns[0] != "time":
            time_column = result.pop("time")
            result.insert(0, "time", time_column)
        return result
    except Exception as exc:
        raise SimulationError(f"Bioscrape simulation failed: {exc}") from exc
    finally:
        sbml_path.unlink(missing_ok=True)

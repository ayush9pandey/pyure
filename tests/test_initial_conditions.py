import pandas as pd

from pyure.initial_conditions import (
    dataframe_to_initial_conditions,
    default_initial_conditions_for_species,
)


def test_default_initial_conditions_filters_to_model_species():
    defaults = default_initial_conditions_for_species(["ATP", "DNA", "not_real"])
    assert defaults["ATP"] == 3750.0
    assert defaults["DNA"] == 0.005
    assert "not_real" not in defaults


def test_dataframe_to_initial_conditions_omits_zero_values():
    dataframe = pd.DataFrame(
        [
            {"Species": "ATP", "Initial concentration": 10.0},
            {"Species": "GTP", "Initial concentration": 0.0},
        ]
    )
    assert dataframe_to_initial_conditions(dataframe) == {"ATP": 10.0}

import pytest

from pyure.exceptions import ModelNotImplementedError
from pyure.model_builders import ModelHierarchy, MostDetailedBuilder, generate_model


def test_modular_hierarchy_generates_a_crn():
    generated = generate_model(ModelHierarchy.PROTEIN_ONLY, "ATGCAA")
    assert generated.species_names
    assert "protein_GFP" in generated.species_names
    assert generated.reaction_count > 0
    assert generated.sbml_xml.startswith("<?xml")


def test_lumped_resources_uses_basic_pure_shape():
    generated = generate_model(ModelHierarchy.LUMPED_ENERGY_RESOURCES, "ATGCAA")
    assert "metabolite_NTPs" in generated.species_names
    assert "metabolite_AAs" in generated.species_names
    assert generated.summary["Builder"] == "BasicPURE"


def test_phenomenological_is_skipped_for_now():
    with pytest.raises(ModelNotImplementedError, match="skipped"):
        generate_model(ModelHierarchy.PHENOMENOLOGICAL, "ATGCAA")


def test_most_detailed_builder_injects_requested_sequence():
    script = "dna_seq= 'ATGAAA'\nprint(dna_seq)\n"
    updated = MostDetailedBuilder._inject_sequence(script, "ATGCAA")
    assert "dna_seq = 'ATGCAA'" in updated

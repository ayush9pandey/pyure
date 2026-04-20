"""Model hierarchy registry and generation wrappers."""

from __future__ import annotations

import contextlib
import csv
import io
import os
import re
import time
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Protocol

from pyure.exceptions import ModelGenerationError, ModelNotImplementedError
from pyure.export import crn_to_sbml_xml
from pyure.initial_conditions import default_initial_conditions_for_species
from pyure.validation import coding_region_summary, validate_dna_sequence


class ModelHierarchy(str, Enum):
    PROTEIN_ONLY = "protein_only"
    TX_TL_ONLY = "tx_tl_only"
    PHENOMENOLOGICAL = "phenomenological"
    LUMPED_ENERGY_RESOURCES = "lumped_energy_resources"
    NUCLEOTIDE = "nucleotide"


HIERARCHY_OPTIONS: list[tuple[ModelHierarchy, str]] = [
    (ModelHierarchy.PROTEIN_ONLY, "protein expression only"),
    (ModelHierarchy.TX_TL_ONLY, "TX-TL only model"),
    # Phenomenological level intentionally hidden for now.
    # (ModelHierarchy.PHENOMENOLOGICAL, "phenomenological"),
    (
        ModelHierarchy.LUMPED_ENERGY_RESOURCES,
        "lumped resources",
    ),
    (ModelHierarchy.NUCLEOTIDE, "nucleotide-level detail"),
]


@dataclass(frozen=True)
class GeneratedModel:
    hierarchy: ModelHierarchy
    label: str
    dna_sequence: str
    crn: object
    sbml_xml: str
    species_names: list[str]
    reaction_count: int
    build_seconds: float
    summary: dict[str, object]
    initial_conditions: dict[str, float]
    build_log: str = ""


class ModelBuilder(Protocol):
    hierarchy: ModelHierarchy
    label: str

    def build(self, dna_sequence: str) -> GeneratedModel:
        """Build a model for the given sequence."""


class MostDetailedBuilder:
    hierarchy = ModelHierarchy.NUCLEOTIDE
    label = "nucleotide-level detail"

    def __init__(self, script_path: Path | None = None):
        repo_root = Path(__file__).resolve().parents[2]
        self.script_path = script_path or repo_root / "pure_most_detailed.py"

    def build(self, dna_sequence: str) -> GeneratedModel:
        sequence = validate_dna_sequence(dna_sequence)
        if not self.script_path.exists():
            raise ModelGenerationError(
                f"Detailed model script was not found at {self.script_path}."
            )

        script_text = self.script_path.read_text(encoding="utf-8")
        script_text = self._inject_sequence(script_text, sequence)
        namespace: dict[str, object] = {
            "__name__": "_pyure_pure_most_detailed_runtime",
            "__file__": str(self.script_path),
            "csv": csv,
        }

        output = io.StringIO()
        start = time.perf_counter()
        original_cwd = Path.cwd()
        try:
            os.chdir(self.script_path.parent)
            with contextlib.redirect_stdout(output):
                exec(compile(script_text, str(self.script_path), "exec"), namespace)
        except Exception as exc:
            raise ModelGenerationError(
                f"Failed to generate the detailed PURE CRN: {exc}"
            ) from exc
        finally:
            os.chdir(original_cwd)

        crn = namespace.get("Combine_PURE")
        if crn is None:
            raise ModelGenerationError(
                "Detailed model script completed but did not create Combine_PURE."
            )

        sbml_xml = crn_to_sbml_xml(crn)
        species_names = _sbml_species_ids(sbml_xml) or _species_names(getattr(crn, "species", []))
        reactions = getattr(crn, "reactions", [])
        build_seconds = time.perf_counter() - start
        coding_summary = coding_region_summary(sequence)
        summary = {
            "DNA length": len(sequence),
            "Coding start index": coding_summary["start_index"],
            "Complete codons from ATG": coding_summary["codons"],
            "Stop codon present": coding_summary["has_stop"],
            "Species": len(species_names),
            "Reactions": len(reactions),
        }
        initial_conditions = default_initial_conditions_for_species(species_names)

        return GeneratedModel(
            hierarchy=self.hierarchy,
            label=self.label,
            dna_sequence=sequence,
            crn=crn,
            sbml_xml=sbml_xml,
            species_names=species_names,
            reaction_count=len(reactions),
            build_seconds=build_seconds,
            summary=summary,
            initial_conditions=initial_conditions,
            build_log=output.getvalue(),
        )

    @staticmethod
    def _inject_sequence(script_text: str, sequence: str) -> str:
        replacement = f"dna_seq = {sequence!r}"
        updated, count = re.subn(
            r"^dna_seq\s*=\s*['\"][ACGTacgt]+['\"]\s*$",
            replacement,
            script_text,
            count=1,
            flags=re.MULTILINE,
        )
        if count != 1:
            raise ModelGenerationError(
                "Could not locate the dna_seq assignment in pure_most_detailed.py."
            )
        return updated


class ModularBioCRNpylerBuilder:
    """Lower-detail BioCRNpyler builder following the pure_mid_detail.py style."""

    def __init__(
        self,
        hierarchy: ModelHierarchy,
        label: str,
        mixture_name: str,
        *,
        use_local_parameters: bool = True,
    ):
        self.hierarchy = hierarchy
        self.label = label
        self.mixture_name = mixture_name
        self.use_local_parameters = use_local_parameters

    def build(self, dna_sequence: str) -> GeneratedModel:
        sequence = validate_dna_sequence(dna_sequence)
        try:
            import biocrnpyler as bcp
        except Exception as exc:
            raise ModelGenerationError(
                "BioCRNpyler is not available in this Python environment."
            ) from exc

        mixture_class = getattr(bcp, self.mixture_name, None)
        if mixture_class is None:
            raise ModelNotImplementedError(
                f"{self.label} needs BioCRNpyler.{self.mixture_name}, "
                "which is not available in this environment."
            )

        coding_summary = coding_region_summary(sequence)
        coding_length = max(3, int(coding_summary.get("codons") or 1) * 3)
        parameters = _modular_parameters(coding_length)
        dna = bcp.DNAassembly(
            name="target",
            promoter="pconst",
            rbs="rbs_strong",
            protein="Protein",
            length=coding_length,
        )

        start = time.perf_counter()
        try:
            mixture_kwargs = {"name": self.hierarchy.value, "components": [dna]}
            if self.use_local_parameters:
                mixture_kwargs["parameters"] = parameters
            mixture = mixture_class(**mixture_kwargs)
            crn = mixture.compile_crn()
        except Exception as exc:
            raise ModelGenerationError(
                f"Failed to generate {self.label} with {self.mixture_name}: {exc}"
            ) from exc

        sbml_xml = crn_to_sbml_xml(crn)
        species_names = _sbml_species_ids(sbml_xml) or _species_names(getattr(crn, "species", []))
        reactions = getattr(crn, "reactions", [])
        build_seconds = time.perf_counter() - start
        initial_conditions = _modular_initial_conditions(species_names)
        summary = {
            "DNA length": len(sequence),
            "Coding start index": coding_summary["start_index"],
            "Complete codons from ATG": coding_summary["codons"],
            "Stop codon present": coding_summary["has_stop"],
            "Species": len(species_names),
            "Reactions": len(reactions),
            "Builder": self.mixture_name,
        }

        return GeneratedModel(
            hierarchy=self.hierarchy,
            label=self.label,
            dna_sequence=sequence,
            crn=crn,
            sbml_xml=sbml_xml,
            species_names=species_names,
            reaction_count=len(reactions),
            build_seconds=build_seconds,
            summary=summary,
            initial_conditions=initial_conditions,
        )


class SkippedBuilder:
    """Hierarchy slot intentionally left out of this pass."""

    def __init__(self, hierarchy: ModelHierarchy, label: str):
        self.hierarchy = hierarchy
        self.label = label

    def build(self, dna_sequence: str) -> GeneratedModel:
        validate_dna_sequence(dna_sequence)
        raise ModelNotImplementedError(
            f"{self.label} is intentionally skipped for now."
        )


def _species_names(species: list[object]) -> list[str]:
    names = []
    for item in species:
        name = getattr(item, "name", None) or str(item)
        names.append(str(name))
    return names


def _sbml_species_ids(sbml_xml: str) -> list[str]:
    """Extract SBML species IDs for Bioscrape-facing names."""
    return re.findall(r"<species\b[^>]*\bid=\"([^\"]+)\"", sbml_xml)


def _modular_parameters(coding_length: int) -> dict[tuple[str | None, str | None, str], float]:
    """Parameter set for BioCRNpyler extract-style hierarchy levels."""
    parameters: dict[tuple[str | None, str | None, str], float] = {}

    def add(mechanism: str, part_ids: list[str | None], values: dict[str, float]) -> None:
        for part_id in part_ids:
            for parameter_name, value in values.items():
                parameters[(mechanism, part_id, parameter_name)] = value

    length = float(max(1, coding_length))
    add("gene_expression", [None, "target", "pconst"], {"kexpress": 0.02})
    add("simple_transcription", [None, "target", "pconst"], {"ktx": 0.05})
    add("simple_translation", [None, "target", "rbs_strong"], {"ktl": 0.1})
    add(
        "transcription_mm",
        [None, "target", "pconst"],
        {"ktx": 0.05, "kb": 100.0, "ku": 10.0},
    )
    add(
        "translation_mm",
        [None, "target", "rbs_strong"],
        {"ktl": 0.1, "kb": 100.0, "ku": 10.0},
    )
    add(
        "rna_degredation_mm",
        [None, "rna_target", "complex_protein_Ribo_rna_target_"],
        {"kdeg": 0.005, "kb": 100.0, "ku": 10.0},
    )
    add("rna_degredation", [None, "rna_target"], {"kdil": 0.005})
    add("one_step_pathway", [None, "target"], {"k": 0.05})
    for mechanism in ("transcription_mm", "translation_mm"):
        add(mechanism, [None, "target", "pconst", "rbs_strong"], {"length": length})
    return parameters


def _modular_initial_conditions(species_names: list[str]) -> dict[str, float]:
    """PURE-informed defaults for modular BioCRNpyler models."""
    available = set(species_names)
    candidates = {
        "dna_target": 0.005,
        "target": 0.005,
        "T7RNAP": 1.0,
        "RNAP": 1.0,
        "RS70S": 3.0,
        "Ribo": 3.0,
        "RNAase": 1.0,
        "ATP": 1.0,
        "NTPs": 1.0,
        "AAs": 1.0,
        "amino_acids": 1.0,
        "Fuel_CP": 1.0,
        "Fuel_3PGA": 1.0,
        "protein_T7RNAP": 1.0,
        "protein_RNAP": 1.0,
        "protein_RS70S": 3.0,
        "protein_Ribo": 3.0,
        "protein_RNAase": 1.0,
        "metabolite_ATP": 1.0,
        "metabolite_NTPs": 1.0,
        "metabolite_AAs": 1.0,
        "metabolite_amino_acids": 1.0,
        "metabolite_Fuel_CP": 1.0,
        "metabolite_Fuel_3PGA": 1.0,
    }
    return {
        species: value
        for species, value in candidates.items()
        if species in available and value != 0
    }


def _basic_pure_class_name() -> str:
    """Return the BioCRNpyler BasicPURE class name used by installed versions."""
    return "BasicPURE"


def get_builder(hierarchy: ModelHierarchy | str) -> ModelBuilder:
    hierarchy_value = ModelHierarchy(hierarchy)
    registry: dict[ModelHierarchy, ModelBuilder] = {
        ModelHierarchy.NUCLEOTIDE: MostDetailedBuilder(),
        ModelHierarchy.LUMPED_ENERGY_RESOURCES: ModularBioCRNpylerBuilder(
            ModelHierarchy.LUMPED_ENERGY_RESOURCES,
            "lumped resources",
            _basic_pure_class_name(),
            use_local_parameters=False,
        ),
        ModelHierarchy.PHENOMENOLOGICAL: SkippedBuilder(
            ModelHierarchy.PHENOMENOLOGICAL,
            "phenomenological",
        ),
        ModelHierarchy.TX_TL_ONLY: ModularBioCRNpylerBuilder(
            ModelHierarchy.TX_TL_ONLY,
            "TX-TL only model",
            "TxTlExtract",
            use_local_parameters=False,
        ),
        ModelHierarchy.PROTEIN_ONLY: ModularBioCRNpylerBuilder(
            ModelHierarchy.PROTEIN_ONLY,
            "protein expression only",
            "ExpressionExtract",
            use_local_parameters=False,
        ),
    }
    return registry[hierarchy_value]


def generate_model(hierarchy: ModelHierarchy | str, dna_sequence: str) -> GeneratedModel:
    """Generate a model for the requested hierarchy."""
    return get_builder(hierarchy).build(dna_sequence)

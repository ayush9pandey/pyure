"""Input validation helpers for PURE CRN model generation."""

from __future__ import annotations

import re

from pyure.exceptions import SequenceValidationError

DNA_ALPHABET = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}


def normalize_dna_sequence(sequence: str) -> str:
    """Return an uppercase DNA sequence with whitespace removed."""
    return re.sub(r"\s+", "", sequence or "").upper()


def validate_dna_sequence(sequence: str) -> str:
    """Validate and normalize a DNA sequence for the current detailed builder."""
    normalized = normalize_dna_sequence(sequence)
    if not normalized:
        raise SequenceValidationError("Enter a DNA sequence before generating a model.")

    invalid = sorted(set(normalized) - DNA_ALPHABET)
    if invalid:
        raise SequenceValidationError(
            "DNA sequence can only contain A, C, G, and T. "
            f"Invalid character(s): {', '.join(invalid)}."
        )

    start = normalized.find("ATG")
    if start == -1:
        raise SequenceValidationError(
            "The detailed translation model needs an ATG start codon."
        )

    coding_region = normalized[start:]
    if len(coding_region) < 6:
        raise SequenceValidationError(
            "The detailed model needs at least one codon after the ATG start codon."
        )

    usable_length = len(coding_region) - (len(coding_region) % 3)
    usable_coding_region = coding_region[:usable_length]
    if len(usable_coding_region) < 6:
        raise SequenceValidationError(
            "The coding region after ATG must include at least two complete codons."
        )

    return normalized


def coding_region_summary(sequence: str) -> dict[str, object]:
    """Return lightweight coding-region details for display."""
    normalized = normalize_dna_sequence(sequence)
    start = normalized.find("ATG")
    if start == -1:
        return {"start_index": None, "codons": 0, "has_stop": False}

    coding_region = normalized[start:]
    usable_length = len(coding_region) - (len(coding_region) % 3)
    codons = [
        coding_region[index : index + 3]
        for index in range(0, usable_length, 3)
    ]
    return {
        "start_index": start,
        "codons": len(codons),
        "has_stop": any(codon in STOP_CODONS for codon in codons),
    }

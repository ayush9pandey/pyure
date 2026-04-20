import pytest

from pyure.exceptions import SequenceValidationError
from pyure.validation import normalize_dna_sequence, validate_dna_sequence


def test_normalize_dna_sequence_removes_whitespace_and_uppercases():
    assert normalize_dna_sequence(" atg\ncaa ") == "ATGCAA"


def test_validate_dna_sequence_rejects_invalid_characters():
    with pytest.raises(SequenceValidationError, match="Invalid character"):
        validate_dna_sequence("ATGBBB")


def test_validate_dna_sequence_requires_start_codon():
    with pytest.raises(SequenceValidationError, match="ATG start codon"):
        validate_dna_sequence("GGGCCCAAATTT")


def test_validate_dna_sequence_accepts_minimal_coding_region():
    assert validate_dna_sequence("GGGATGCAA") == "GGGATGCAA"

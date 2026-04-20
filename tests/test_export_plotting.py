from pyure.export import download_filename
from pyure.plotting import parse_species_input, validate_requested_species


def test_download_filename_sanitizes_slug():
    assert download_filename("nucleotide detail") == "pyure_nucleotide_detail.xml"


def test_parse_species_input_accepts_commas_and_newlines():
    assert parse_species_input("ATP, GTP\nmRNA") == ["ATP", "GTP", "mRNA"]


def test_validate_requested_species_reports_found_and_missing():
    found, missing = validate_requested_species(["ATP", "missing"], ["time", "ATP"])
    assert found == ["ATP"]
    assert missing == ["missing"]

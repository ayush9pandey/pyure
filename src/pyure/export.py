"""SBML export helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path


def crn_to_sbml_xml(crn: object, *, check_validity: bool = True) -> str:
    """Serialize a BioCRNpyler ChemicalReactionNetwork to SBML XML text."""
    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as handle:
        output_path = Path(handle.name)

    try:
        try:
            crn.write_sbml_file(str(output_path), check_validity=check_validity)
        except Exception:
            if not check_validity:
                raise
            crn.write_sbml_file(str(output_path), check_validity=False)
        return output_path.read_text(encoding="utf-8")
    finally:
        output_path.unlink(missing_ok=True)


def download_filename(model_slug: str) -> str:
    """Return a stable SBML download filename for a hierarchy level."""
    safe_slug = "".join(
        char if char.isalnum() or char in {"-", "_"} else "_"
        for char in model_slug.lower()
    ).strip("_")
    return f"pyure_{safe_slug or 'model'}.xml"

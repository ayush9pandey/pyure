"""Rough model-generation time estimates."""

from __future__ import annotations

import math


def estimate_generation_seconds(nucleotide_length: int) -> float:
    """Estimate detailed model generation time from short local benchmarks.

    Local pyure310 smoke timings including wrapper execution and SBML export were
    about 2.6 s at 9 nt and 9.5 s at 48 nt. The long-sequence anchor remains the
    user-observed about 300 s at 800 nt.
    """
    length = max(1, int(nucleotide_length))
    short_length = 9
    mid_length = 48
    long_length = 800
    short_seconds = 2.6
    mid_seconds = 9.5
    long_seconds = 300.0
    if length <= mid_length:
        slope = (mid_seconds - short_seconds) / (mid_length - short_length)
        return max(0.1, short_seconds + slope * (length - short_length))
    exponent = math.log(long_seconds / mid_seconds) / math.log(long_length / mid_length)
    return mid_seconds * (length / mid_length) ** exponent


def format_seconds(seconds: float) -> str:
    """Format seconds into compact UI text."""
    if seconds < 1:
        return f"{seconds:.2f} seconds"
    if seconds < 90:
        return f"{seconds:.1f} seconds"
    return f"{seconds / 60:.1f} minutes"

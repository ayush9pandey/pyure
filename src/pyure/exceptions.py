"""Application-specific exceptions."""


class PyureError(Exception):
    """Base exception for expected application errors."""


class SequenceValidationError(PyureError):
    """Raised when a DNA sequence cannot be modeled."""


class ModelNotImplementedError(PyureError):
    """Raised when a hierarchy level does not have a builder yet."""


class ModelGenerationError(PyureError):
    """Raised when model generation fails."""


class SimulationError(PyureError):
    """Raised when Bioscrape simulation fails."""

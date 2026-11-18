"""seldecomp: raking-based reweighting decomposition of group differences in means."""
from .api import seldecomp
from .errors import SeldecompError, SeldecompWarning
from ._results import SeldecompResults

__version__ = "1.0.0"

__all__ = ["seldecomp", "SeldecompError", "SeldecompResults", "SeldecompWarning", "__version__"]

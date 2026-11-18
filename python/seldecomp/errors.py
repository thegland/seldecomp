"""The exception and the warning category of seldecomp."""
import os
import sys
import warnings

_PACKAGE = os.path.dirname(os.path.abspath(__file__))


class SeldecompError(Exception):
    """An error of seldecomp. rc is the return code the Stata command gives in the same
    situation (198 for an invalid command, 430 when the raking misses its targets, 459 for a
    failure of support, and so on)."""

    def __init__(self, message, rc=198):
        Exception.__init__(self, message)
        self.rc = rc
        self.message = message


class SeldecompWarning(UserWarning):
    """A warning of seldecomp: the text the Stata command prints after "warning: "."""


def warn(message):
    """Issue a SeldecompWarning attributed to the line that called into the package."""
    level = 2
    frame = sys._getframe(1)
    while frame is not None and os.path.dirname(os.path.abspath(frame.f_code.co_filename)) == _PACKAGE:
        frame = frame.f_back
        level += 1
    warnings.warn(message, SeldecompWarning, stacklevel=level)

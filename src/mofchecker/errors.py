# -*- coding: utf-8 -*-
"""Custom error types."""


class NoMetal(KeyError):
    """Error in case there is no metal in structure."""

class Too_many_adding_linkers(KeyError):
    """Error in case there are too many adding linkers."""

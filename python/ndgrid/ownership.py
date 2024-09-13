"""Ownership."""

class Ownership(object):
    """Ownership base class."""


class Owned(Ownership):
    """Owned by the current process."""


class Ghost(Ownership):
    """Owned by the current process."""

    def __init__(self, process: int, index: int):
        """Initialise."""
        self.process = process
        self.index = index

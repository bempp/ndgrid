"""Ownership."""

from abc import ABC, abstractproperty


class Ownership(ABC):
    """Ownership base class."""

    @abstractproperty
    def is_owned(self) -> bool:
        """Is this entity owned by the current process?"""


class Owned(Ownership):
    """Owned by the current process."""

    @property
    def is_owned(self) -> bool:
        """Is this entity owned by the current process?"""
        return True


class Ghost(Ownership):
    """Owned by the current process."""

    def __init__(self, process: int, index: int):
        """Initialise."""
        self.process = process
        self.index = index

    @property
    def is_owned(self) -> bool:
        """Is this entity owned by the current process?"""
        return False

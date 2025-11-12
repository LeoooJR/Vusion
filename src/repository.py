"""
A module to manage the repositories of objects.
"""

from abc import ABC, abstractmethod
from typing import Final, Generic, Iterator, TypeVar

T: TypeVar = TypeVar("T")  # Type of the objects in the repository


class Repository(ABC, Generic[T]):
    """
    A repository of objects.
    """

    def __init__(self) -> None:
        """Initialize the repository"""

        # Dictionary to store the objects in the repository
        self.repository: Final[dict[str, T]] = {}

    @abstractmethod
    def populate(self, items: list[T]) -> None:
        """
        Populate the repository with a list of objects.

        Args:
            items (list[T]): The list of objects to add.
        """

    @abstractmethod
    def add(self, item: T) -> None:
        """
        Add an object to the repository.

        Args:
            item (T): The object to add.
        """

    @abstractmethod
    def remove(self, item: T) -> None:
        """
        Remove an object from the repository.

        Args:
            item (T): The object to remove.
        """

    @abstractmethod
    def __iter__(self) -> Iterator[tuple[str, T]]:
        """
        Iterate over the repository.

        Yields:
            tuple[str, T]: A tuple containing the key and the object.
        """

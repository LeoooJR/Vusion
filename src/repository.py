from abc import ABC, abstractmethod
from typing import Generic, Iterator, TypeVar

T = TypeVar("T")


class Repository(ABC, Generic[T]):
    """
    A repository of objects.
    """

    @abstractmethod
    def populate(self, items: list[T]):
        """
        Populate the repository with a list of objects.
        """

    @abstractmethod
    def add(self, item: T):
        """
        Add an object to the repository.
        """

    @abstractmethod
    def remove(self, item: T):
        """
        Remove an object from the repository.
        """

    @abstractmethod
    def __iter__(self) -> Iterator[T]:
        """
        Iterate over the repository.
        """

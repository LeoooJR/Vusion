from abc import ABC, abstractmethod
from typing import Generic, TypeVar

T = TypeVar('T')

class Repository(ABC, Generic[T]):

    @abstractmethod
    def populate(self, items: list[T]):

        pass

    @abstractmethod
    def add(self, item: T):

        pass

    @abstractmethod
    def remove(self, item: T):

        pass

    @abstractmethod
    def __iter__(self):

        pass
from abc import ABC, abstractmethod

class BaseFilter(ABC):
    def __init__(self, order, fc):
        self.order = order
        self.fc = fc  # frecuencia de corte
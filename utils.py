from dataclasses import dataclass

@dataclass
class Pars:
    N: int
    rh: float
    alpha: float
    beta: float
    eta: float
    nx: int
    ny: int
    nz: int
    boxsize: float
    repeat: int
    checkerror: bool

    def __iter__(self):
        # Iterate over field names and values
        for field in self.__annotations__:
            yield field, getattr(self, field)
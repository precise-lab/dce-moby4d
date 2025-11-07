from enum import Enum

class Chromophore(Enum):
    HB = 0
    HB02 = 1
    WATER = 2
    LIPID = 3
    MELATONIN = 4
    CA    = 5
    NCOMP = 6

class ChromophoreDecomposition:
    def __init__(self, spectra):
        self.spectra = spectra
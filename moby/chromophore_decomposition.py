import os
from enum import Enum
import h5py
import numpy as np

class Chromophore(Enum):
    HB = 0
    HBO2 = 1
    WATER = 2
    LIPID = 3
    MELATONIN = 4
    CA    = 5

class ChromophoreDecomposition:
    def __init__(self, wavelength, spectra, units):
        self.wavelength = wavelength
        self.spectra = spectra
        self.units    = units
        self.c_thb_b  = None

    def setThBConcentration(self, c_thb_b: float):
        self.c_thb_b  = c_thb_b #M

    @classmethod
    def fromFile(cls, chromophores_var: dict = None, fname: str=None, *, verbose=False):
        if fname == None:
            fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spectra.h5")
        if chromophores_var == None:
            chromophores_var = {}
            chromophores_var[Chromophore.HB] = "e_hb"
            chromophores_var[Chromophore.HBO2] = "e_hbO2"
            chromophores_var[Chromophore.WATER] = "mu_a_water"
            chromophores_var[Chromophore.CA]   = "e_icg_6.5e-6_M"


        with h5py.File(fname, 'r') as fid:
            if verbose:
                print("Reading spectra from: ", fname)
            colnames = fid["colnames"]
            all_units    = fid["units"]
            data  = fid["spectra"]
            colnames = [c.decode("utf-8") for c in colnames]
            all_units = [u.decode("utf-8") for u in all_units]
            if verbose:
                for i in range(len(colnames) ):
                    print(colnames[i], " ", all_units[i])

            assert(colnames[0] == 'wavelength')
            wavelength = data[:,0]
            spectra = {}
            units   = {}

            for c in chromophores_var:
                idx = colnames.index( chromophores_var[c] )
                spectra[c] = data[:,idx]
                units[c] = all_units[idx]

            return cls(wavelength, spectra, units)
        
    def get_pure_mu_a(self, wavelength: float, c_list: list[Chromophore]):
        out = {}
        for c in c_list:
            s = np.interp(wavelength, self.wavelength, self.spectra[c])
            if c in [Chromophore.HB, Chromophore.HBO2]:
                out[c] = np.log(10.)*self.c_thb_b*s
            elif c == Chromophore.CA:
                out[c] = np.log(10.)*s
            else:
                out[c] = s

        return out



            


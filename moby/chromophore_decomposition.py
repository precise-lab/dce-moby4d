import os
from enum import Enum
import h5py
import numpy as np

class Chromophore(Enum):
    HB = 0
    HB02 = 1
    WATER = 2
    LIPID = 3
    MELATONIN = 4
    CA    = 5
    NC    = 6

class ChromophoreDecomposition:
    def __init__(self, wavelength, spectra, units):
        self.wavelength = wavelength
        self.spectra = spectra
        self.units    = units

    @classmethod
    def fromFile(cls, chromophores_var: dict, fname: str=None):
        if fname == None:
            fname = os.path.join(os.path.abspath(__file__), "spectra.h5")

        with h5py.File(fname, 'r') as fid:
            print("Reading spectra from: ", fname)
            colnames = fid["colnames"]
            all_units    = fid["units"]
            data  = fid["spectra"]
            for i in range(colnames.shape[0]):
                print(colnames[i], " ", all_units[i])

            assert(colnames[0] == 'wavelength')
            wavelength = data[:,0]
            spectra = np.zeros((wavelength.shape[0], Chromophore.NC), dtype = np.float64)
            units   = np.zeros(Chromophore, dtype=all_units.dtype)

            for c in chromophores_var:
                idx = (colnames == chromophores_var[c] )
                assert np.any(idx)
                spectra[:,c] = data[:,idx]
                units[c] = all_units[idx]

            return cls(wavelength, spectra, units)
        
    def get_pure_mu_a(self, wavelength: float, c_thb_b: float, c_CA_b: float, c_list: list[Chromophore]):
        out = np.array(len(c_list))
        s = self.spectra[ self.wavelength==wavelength ] #FIXME use interpolation
        for c, i in enumerate(c_list):
            if c in [Chromophore.HB, Chromophore.HB02]:
                out[i] = np.log(10.)*c_thb_b*s[c]
            elif c == Chromophore.CA:
                out[i] = np.log(10.)*c_CA_b*s[c]
            else:
                out[i] = s[c]

        return out



            


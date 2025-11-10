import numpy as np
import os
import csv




class TissueComposition:
    def __init__(self,
                 volume_fractions: dict,
                 ref_reduced_scattering: dict,
                 ref_wavelength: dict,
                 b_exponent: dict,
                 tissue2label: dict,
                 label2glabel: dict):
        
        self.volume_fractions   = volume_fractions
        self.ref_reduced_scattering = ref_reduced_scattering
        self.ref_wavelength     = ref_wavelength
        self.b_exponent         = b_exponent
        self.tissue2label       = tissue2label
        self.label2glabel       = label2glabel

        # Normal hemoglobin mass concentration ranges of mice
        self.mc_hb = 150 # [g/L]

        # Molecular weight of mammalian hemoglobin
        self.mw_hb = 64500 # [g/mol]


    @property
    def c_thb_b(self):
        "Hemoglobin concentration of blood"
        return self.mc_hb/self.mw_hb # [mol/L, M]
    
    @classmethod
    def create(cls):
        fpath=os.path.abspath(__file__)
        fdir = os.path.dirname(fpath)
        fname = os.path.join(fdir, 'tissue_composition.csv')
        return cls.createFromFile(fname)
        
    @classmethod
    def createFromFile(cls, fname):
        fid = open(fname, newline='')
        reader = csv.DictReader(fid)
        
        volume_fractions = {}
        ref_reduced_scattering = {}
        ref_wavelength = {}
        b_exponent  = {}
        tissue2label = {}
        label2glabel = {}
        for row in reader:
            tissue = row['tissue']
            label = int(row['label'])
            group_label = int(row['group_label'])
            tissue2label[tissue] = label
            label2glabel[label] = group_label
            volume_fractions[label] = {"blood": float( row['blood'] ), "water": float( row['water'] )}
            ref_reduced_scattering[label] = float( row['mu_s_ref'] ) * (1. - float(row['g']))
            ref_wavelength[label] = float(row['ref_wavelength'])
            b_exponent = float(row['b_exponent'])

        return cls(volume_fractions, ref_reduced_scattering, ref_wavelength, b_exponent, tissue2label, label2glabel)
    
    def volumeFractionMap(self, label_map: np.array, comp: str):

        unique_labels = np.unique(label_map[:])
        for ul in unique_labels[:]:
            assert ul in self.volume_fractions

        vf = np.zeros_like(label_map, dtype=np.float64)
        for (label, val) in self.volume_fractions.items():
            vf[label_map==label] = val[comp]

        return vf
    
    def getVolumeFraction(self, label: int, comp: str):
        return self.volume_fractions[label][comp]
    
    def reducedScatteringMap(self, label_map: np.array, wavelength: float):
        unique_labels = np.unique(label_map[:])
        for ul in unique_labels[:]:
            assert ul in self.ref_reduced_scattering

        rs = np.zeros_like(label_map, dtype=np.float64)
        for label in self.ref_reduced_scattering.keys():
            rs[label_map==label] = self.getReducedScattering(label, wavelength)

    def getReducedScattering(self, label: int, wavelength: float):
            rrs = self.ref_reduced_scattering[label]
            rw  = self.ref_wavelength[label]
            b = self.b_exponent[label]
            return rrs*( (rw/wavelength)**b )

    def gLabelMap(self, label_map: np.array):

        unique_labels = np.unique(label_map[:])
        for ul in unique_labels[:]:
            assert ul in self.label2glabel

        glabel_map = np.zeros_like(label_map)
        for (label, glabel) in self.label2glabel.items():
            glabel_map[label_map==label] = glabel

        return glabel_map


import numpy as np
import os
import csv
from enum import Enum

class Component(Enum):
    BLOOD = 0
    WATER = 1
    NCOMP = 2

class TissueComposition:
    def __init__(self, volume_fractions: dict, tissue2label: dict, label2glabel: dict):
        self.volume_fractions = volume_fractions
        self.tissue2label     = tissue2label
        self.label2glabel     = label2glabel

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
        tissue2label = {}
        label2glabel = {}
        for row in reader:
            tissue = row['tissue']
            label = int(row['label'])
            group_label = int(row['group_label'])
            tissue2label[tissue] = label
            label2glabel[label] = group_label
            volume_fractions[label] = (row['blood'], row['water'])

        return cls(volume_fractions, tissue2label, label2glabel)
    
    def volumeFractionMap(self, label_map: np.array, comp: Component):

        unique_labels = np.unique(label_map[:])
        for ul in unique_labels[:]:
            assert ul in self.volume_fractions

        vf = np.zeros_like(label_map, dtype=np.float64)
        for (label, val) in self.volume_fractions.items():
            vf[label_map==label] = val[comp]

        return val

    def gLabelMap(self, label_map: np.array):

        unique_labels = np.unique(label_map[:])
        for ul in unique_labels[:]:
            assert ul in self.label2glabel

        glabel_map = np.zeros_like(label_map)
        for (label, glabel) in self.label2glabel.items():
            glabel_map[label_map==label] = glabel

        return glabel_map


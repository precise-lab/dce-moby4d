import dolfin as dl
from .tissue_composition import TissueComposition
from .pkModel import PKModel
from .chromophore_decomposition import ChromophoreDecomposition

class FEMPhantom:
    def __init__(self, tissue_composition: TissueComposition,
                       pk_model: PKModel)
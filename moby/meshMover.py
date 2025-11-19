import dolfin as dl
import numpy as np
import h5py

class MeshMover():
    def __init__(self, mesh, fname):
        self.mesh = mesh
        self.orig_coords = np.copy(self.mesh.coordinates())
        self.fname = fname

    def move(self, time):
        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            dt = dset.attrs["dt"]
            h = dset.attrs["spacing"]


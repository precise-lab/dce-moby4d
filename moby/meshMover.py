import dolfin as dl
import ufl
import numpy as np
import h5py

import hippylib as hp

class MeshMover():
    def __init__(self, mesh, fname, verbose=False):
        self.mesh = mesh
        self.orig_coords = np.copy(self.mesh.coordinates())
        self.fname = fname
        self.ref_Vh = dl.VectorFunctionSpace(self.mesh, "CG", 1)
        self.verbose = verbose

    def dt(self):
        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            dt = dset.attrs["dt"]
        
        return dt
        
    def nframes(self):

        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            Nx, Ny, Nz, ncomp, nframes = dset.shape

        return nframes

        

    def move(self, time):

        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            dt = dset.attrs["dt"]
            h = dset.attrs["spacing"]

            Nx, Ny, Nz, ncomp, nframes = dset.shape
            assert ncomp == 3
            iframe = int( np.floor(time/dt) ) % nframes

            disp_array = dset[:,:,:,:, iframe]

            
        disp_exp = hp.NumpyVectorExpression3D(3)
        disp_exp.setData(disp_array, *h[:])
        d = dl.interpolate(disp_exp, self.ref_Vh)

        if self.verbose:
            d_np = np.abs( d.vector().gather_on_zero() )
            d_np = np.reshape(d_np, (int(d_np.shape[0]/3), 3))
            nd_np = np.sqrt(np.sum(d_np*d_np, axis=1))
            if nd_np.size > 0:
                print("-->",  np.max(nd_np), " ", np.max(d_np[:,0]), " ", np.max(d_np[:,1]), " ", np.max(d_np[:,2]))

        self.mesh.coordinates()[:] = self.orig_coords
        dl.ALE.move(self.mesh, d)

        DG0 = dl.FunctionSpace(self.mesh, "DG", 0)
        test = dl.TestFunction(DG0)
        vol = dl.assemble(test*ufl.dx)
        vol_np = vol.gather_on_zero()
        if vol_np.size > 0:
            print( np.min(vol_np), " ", np.max(vol_np) )




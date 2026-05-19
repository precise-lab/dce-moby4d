import dolfin as dl
import ufl
import numpy as np
import scipy.ndimage
import h5py

import hippylib as hp

class MeshMover():
    def __init__(self, mesh, fname, verbose=False):
        self.mesh = mesh
        self.orig_coords = np.copy(self.mesh.coordinates())
        self.fname = fname
        self.ref_Vh = dl.VectorFunctionSpace(self.mesh, "CG", 1)
        self.vDG0 = dl.VectorFunctionSpace(self.mesh, "DG", 0)
        self.ref_Vh_scalar = dl.FunctionSpace(self.mesh, "DG", 0)
        self.verbose = verbose

        self._dt = None
        self._spacing = None
        self._nframes = None

    def _fetch_attrs(self):
        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            self._dt = dset.attrs["dt"]
            self._spacing = dset.attrs["spacing"]
            _, _, _, _, self._nframes = dset.shape

    @property
    def dt(self):
        if self._dt is None:
            self._fetch_attrs()

        return self._dt

    
    @property
    def nframes(self):
        if self._nframes is  None:
            self._fetch_attrs()

        return self._nframes
    
    @property
    def spacing(self):
        if self._spacing is None:
            self._fetch_attrs()

        return self._spacing

    def get_imposed(self,time):

        dt = self.dt
        h = self.spacing
        nframes = self._nframes
        iframe = int( np.floor(time/dt) ) % nframes

        with h5py.File(self.fname) as fid:
            iset = fid["interpolated"]
            interp_array = iset[:,:,:, iframe]

        imposed_array = (interp_array==0)
        imposed_array = imposed_array.astype(np.float64)
        imposed_exp =  hp.NumpyScalarExpression3D()
        imposed_exp.setData(imposed_array, *h[:])
        imposed = dl.interpolate(imposed_exp, self.ref_Vh_scalar)
        imposed.rename("Imposed", "Imposed")

        return imposed

    def compute_displacement2(self, time, dx):

        disp_DG0 = self.get_displacement(time, self.vDG0, n_smoothing_steps=3)

        # lung, airways, liver, gal-bladder, spleen, kidney, st_wall, st_cnts, spine, rib
        ids = [22, 23, 29, 30, 37, 34, 31, 32, 82]
        # spine, femor
        ids_still = [94, 88]

        uh, vh = dl.TrialFunction(self.ref_Vh), dl.TestFunction(self.ref_Vh)

        Aform = ( ufl.inner( ufl.sym(ufl.grad(uh)), ufl.sym(ufl.grad(vh)))*dx 
                 + ufl.inner(ufl.div(uh), ufl.div(vh))*dx )
        
        disp = dl.Function(self.ref_Vh, name="Displacement")
        bform = ufl.inner(disp, vh)*dl.Constant(10)*dx
                 
        for id in ids:
            Aform += ufl.inner(uh, vh)*dl.Constant(10)*dx(id)
            bform += ufl.inner(disp_DG0, vh)*dl.Constant(10)*dx(id)

        for id in ids_still:
            Aform += ufl.inner(uh, vh)*dl.Constant(10)*dx(id)

        A, b = dl.assemble_system(Aform, bform)
        solver = dl.PETScKrylovSolver("cg", "petsc_amg")

        solver.solve(A, disp.vector(), b)

        return disp

    def compute_displacement(self, time):
        imposed = self.get_imposed(time)
        disp_DG0 = self.get_displacement(time, self.vDG0, n_smoothing_steps=3)

        uh, vh = dl.TrialFunction(self.ref_Vh), dl.TestFunction(self.ref_Vh)

        Aform = ( ufl.inner( ufl.sym(ufl.grad(uh)), ufl.sym(ufl.grad(vh)))*ufl.dx 
                 + ufl.inner(ufl.div(uh), ufl.div(vh))*ufl.dx
                 + ufl.inner(uh, vh)*dl.Constant(10)*imposed*ufl.dx )
        
        bform = ufl.inner(disp_DG0, vh)*dl.Constant(10)*imposed*ufl.dx

        disp = dl.Function(self.ref_Vh, name="Displacement")

        A, b = dl.assemble_system(Aform, bform)
        solver = dl.PETScKrylovSolver("cg", "petsc_amg")

        solver.solve(A, disp.vector(), b)

        return disp

    def get_disp_array(self, time, n_smooth=0):
        with h5py.File(self.fname) as fid:
            dset = fid["motion"]
            dt = dset.attrs["dt"]

            _, _, _, ncomp, nframes = dset.shape
            assert ncomp == 3
            iframe = int( np.floor(time/dt) ) % nframes

            disp_array = dset[:,:,:,:, iframe]
            iset = fid["interpolated"]
            interp_array = iset[:,:,:, iframe]

        imposed_array = (interp_array==0)
        requested_shape = interp_array.shape+(ncomp,)
        
        imposed_array = np.broadcast_to(imposed_array[:,:,:,np.newaxis], requested_shape  )

        for _ in range(n_smooth):
            disp_array_new = scipy.ndimage.uniform_filter(disp_array, size=[3,3,3,1])
            disp_array_new[imposed_array] = disp_array[imposed_array]
            disp_array = disp_array_new

        return disp_array

    def get_displacement(self, time, Vh=None, n_smoothing_steps=1): 

        if Vh is None:
            Vh = self.ref_Vh

        disp_array = self.get_disp_array(time, n_smoothing_steps)
        h = self.spacing

        disp_exp = hp.NumpyVectorExpression3D(3)
        disp_exp.setData(disp_array, *h[:])
        d = dl.interpolate(disp_exp, Vh)
        d.rename("Displacement", "Displacement")

        return d

    def move(self, time):

        d = self.compute_displacement(time)

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

    def move2(self, time, dx):

        d = self.compute_displacement2(time, dx)

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




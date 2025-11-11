import dolfin as dl
import numpy as np
import hippylib as hp

class CartesianGridResampler:
    def __init__(self, Vh: dl.FunctionSpace, fov: np.array, N: np.array, subsampling: int):
        #fov = [xmin, xmax; ymin, ymax; zmin, zmax]
        # N = [Nx, Ny, Nz]
        self.Vh = Vh
        self.fov = fov
        self.N = N
        self.h = (fov[:,1] - fov[:,0])/self.N

        if subsampling == 0:
            offset = np.zeros(3)
            self.B = [self._getB(offset)]
        elif subsampling == 1:
            oxyz = [self.h[i]*np.array([-0.5/np.sqrt(3.), 0.5/np.sqrt(3.)]) for i in range(3)] 
            OX, OY, OZ = np.meshgrid(*oxyz)
            offsets = np.hstack([np.reshape(OXYZ,  (OXYZ.size,1)) for OXYZ in [OX, OY, OZ]])

            self.B = [self._getB(offsets[i,:]) for i in range(1,offsets.shape[0])]
 

    def _getB(self,  offset):
        xyz = [self.fov[i,0]+offset[i] + np.linspace(.5, (self.N[i]-.5), self.N[i])*self.h[i] for i in range(3)]
        XX, YY, ZZ = np.meshgrid(*xyz)
        points = np.hstack([np.reshape(XYZ, (XYZ.size,1)) for XYZ in [XX,YY,ZZ] ])

        return hp.assemblePointwiseObservation(self.Vh, points)


    def __call__(self, f: dl.Function):
        rescale = 1./len(self.B)
        f_np = rescale*(self.B[0]*f.vector()).gather_on_zero()
        for B in self.B[1:]:
            f_np = f_np + rescale*(B*f.vector()).gather_on_zero()
            
        if f_np.size > 0:
            return np.reshape(f_np, self.N).astype(np.float32)
        else:
            return None
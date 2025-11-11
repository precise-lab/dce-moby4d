import dolfin as dl
import ufl
import numpy as np
import scipy.io as io

import pyvista as pv
import matplotlib.pyplot as plt

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp



sys.path.append("../")
import moby
import time



def move(mesh):
    Vh = dl.VectorFunctionSpace(mesh, "CG", 1)
    d  = dl.interpolate(dl.Expression(("0.", "0.", "-2.*x[0]*(1-x[0])*x[1]*(1.-x[1])"), degree=1), Vh)
    dl.ALE.move(mesh, d)

if __name__=="__main__":
    comm = dl.MPI.comm_world
    N = [200, 200, 100]
    mesh = dl.UnitCubeMesh(*N)
    move(mesh)
    Vh = dl.FunctionSpace(mesh, "DG", 0)
    f = dl.interpolate( dl.Expression("-1.", degree=1), Vh)

    fov = np.array([[0., 1.]]*3)
    resampler = moby.CartesianGridResampler(Vh, fov, N, 0)
    f_np = resampler(f)

    with dl.XDMFFile("tcgr.xdmf") as fid:
        fid.write(f)



    if 0 == comm.rank:
        plt.imshow(f_np[:,:,-3])
        plt.savefig('tcgr2d.png')

        pl = pv.Plotter()
        pl.add_volume(f_np, cmap='viridis')
        pl.show(screenshot='tcgr.png')
    

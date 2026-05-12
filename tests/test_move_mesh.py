import dolfin as dl
import ufl
import numpy as np
import pyvista as pv
import h5py
import matplotlib.pyplot as plt

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp

sys.path.append("../")
import moby

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='test_move_mesh', fromfile_prefix_chars='@')
    parser.add_argument('-m', '--mesh', default = "moby_mesh.xdmf")
    parser.add_argument('-d', '--disp', default = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors.h5')
    
    args = parser.parse_args()

    comm = dl.MPI.comm_world

    imaging_time = 0.
    fov = np.array([[0., 37.2], [0., 37.2], [26.25, 56.25]])
    N = [248, 248, 200]

    mesh = dl.Mesh(comm)

    with dl.XDMFFile(args.mesh) as fid:
        fid.read(mesh)
        geo_dim = mesh.geometry().dim()
        c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
        fid.read(c_labels, "c_labels")

    dx   = ufl.Measure("dx", domain = mesh, subdomain_data=c_labels, metadata={"quadrature_degree": 0})

    mover = moby.MeshMover(mesh, args.disp, verbose=True)

    dt = mover.dt()
    nframes = mover.nframes()
    times = np.arange(0, dt*nframes - 0.5*dt, dt)
    pos = np.zeros((nframes, 3))

    for i in np.arange(nframes):
         t = i*dt
         times[i] = t
         if comm.rank == 0:
             print(f"Time {t} s")
         mover.move(t)

         x = ufl.SpatialCoordinate(mesh)
         vol = dl.assemble(dl.Constant(1.)*dx(97)+dl.Constant(1.)*dx(98))
         for jcoor in range(3):
            pos[i,jcoor] = dl.assemble(x[jcoor]*dx(97) + x[jcoor]*dx(98))/vol

    if comm.rank == 0:
        plt.plot(times, pos[:,0]-pos[0,0], label = 'X-pos')
        plt.plot(times, pos[:,1]-pos[0,1], label = 'Y-pos')
        plt.plot(times, pos[:,2]-pos[0,2], label = 'Z-pos')
        plt.legend()
        plt.show()


 
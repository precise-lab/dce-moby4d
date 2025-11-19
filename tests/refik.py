"""
Quick code for Refik to generate numpy maps of mu_a and mu_ps at three different time steps and wavelengths

First download from Box: /Virtual_Imaging_Trials/Mouse/MOBY_phantom_motion/phantom_anatomy_0.h5 and save it in the current folder

Next generate a fitted structured mesh with
> python generate_mesh_structured.py --zmin 100 --zmax 450 --output moby_mesh_100-450.xdmf

Then run
> mpirun -n 8 python refik.py

It output a file called muamusp.h5 with two fields:
- mu_sp of size Nx x Ny x Nz x 1 x number_of_wavelengths
- mu_a  of size Nx x Ny x Nz x number_of_time_snapshots x number_of_wavelengths

where
 - Nx, Ny, Nz denote the number of voxels in the x, y, and z directions;
 - number_of_time_snapshots is the number of time snapshots considered (3: before injection, AIF peak, wash-out)
 - number_of_wavelengths is the number of used wavelengths 8
"""

import dolfin as dl
import ufl
import numpy as np
import pyvista as pv
import h5py

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp

sys.path.append("../")
import moby

def get_boundingbox(Vh):
      x = dl.interpolate(dl.Expression("x[0]", degree=1), Vh)
      y = dl.interpolate(dl.Expression("x[1]", degree=1), Vh)
      z = dl.interpolate(dl.Expression("x[2]", degree=1), Vh)

      bb = np.array([ [cc.vector().min(), cc.vector().max()] for cc in [x,y,z] ] )

      return bb

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='test_fluence', fromfile_prefix_chars='@')
    parser.add_argument('-m', '--mesh', default = "moby_mesh_100-450.xdmf")
    
    args = parser.parse_args()

    comm = dl.MPI.comm_world

    wavelengths = [730, 750, 770, 780, 790, 800, 820, 840] #nm
    imaging_times = [-30., 43.2, 600] #seconds
    N = [248, 248, 350]

    tissueComposition = moby.TissueComposition.create()
    chromophores = [moby.Chromophore.HB, moby.Chromophore.HBO2, moby.Chromophore.WATER, moby.Chromophore.CA]

    time = np.arange(-30., 450.01, 0.05)
    Ktrans  = 0.3/60 # 1/seconds
    kep     = 0.75/60  # 1/seconds  
    aif = moby.AIF()
    pkModel = moby.PKModel(aif, time, Ktrans, kep)
    chromophoresDec = moby.ChromophoreDecomposition.fromFile(verbose=(0==comm.rank))
    chromophoresDec.setThBConcentration(tissueComposition.c_thb_b)


    mesh = dl.Mesh(comm)

    with dl.XDMFFile(args.mesh) as fid:
        fid.read(mesh)
        geo_dim = mesh.geometry().dim()
        c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
        fid.read(c_labels, "c_labels")

    Vh_phi = dl.FunctionSpace(mesh, "CG", 1)
    Vh_m  = dl.FunctionSpace(mesh, "DG", 0)

    fov = get_boundingbox(Vh_phi)
    if 0==comm.rank:
        print(fov)

    resampler = moby.CartesianGridResampler(Vh_m, fov, N, subsampling=1)
        
    femPhantom = moby.FEMPhantom(Vh_m, Vh_phi, c_labels, tissueComposition, chromophoresDec, pkModel)
    sat_map = {"artery": 0.98, "vein": 0.7, "tumor": .371, "tumor_core": 0}
    so2 = femPhantom.compute_oxygen_saturation(sat_map )

    mu_sp_np = np.zeros(N+[1,len(wavelengths)])
    mu_a_np = np.zeros(N+[len(imaging_times),len(wavelengths)])

    for i, wavelength in enumerate(wavelengths):
        mu_sp = femPhantom.compute_mu_sp(wavelength)
        mu_sp_np[:,:,:,0,i] = resampler(mu_sp)
        for j, imaging_time in enumerate(imaging_times):
            mu_a = femPhantom.compute_mu_a(wavelength, imaging_time, so2)
            mu_a_np[:,:,:, j, i] = resampler(mu_a) 

    if 0 == comm.rank:
        pl = pv.Plotter(shape=(1, 3), border=False)
        for j, imaging_time in enumerate(imaging_times):
            pl.subplot(0,j)
            pl.add_volume(np.squeeze(mu_a_np[:,:,:,j, wavelengths.index(800)]), cmap='viridis', opacity="sigmoid")
            pl.add_title("Time = {}s".format(imaging_time))
        pl.show(screenshot='mua.png')

    h = (fov[:,1] - fov[:,0])/np.array(N)
    if 0 == comm.rank:
        with h5py.File("muamusp.h5", "w") as fid:
             mu_sp_set = fid.create_dataset("mu_sp", data = mu_sp_np, compression="lzf")
             mu_sp_set.attrs["spacing"] = h
             mu_sp_set.attrs["spacing_units"] = "mm"
             mu_sp_set.attrs["fov"] = fov
             mu_sp_set.attrs["wavelength"] = np.array(wavelengths)
             mu_a_set = fid.create_dataset("mu_a", data = mu_a_np, compression="lzf")
             mu_a_set.attrs["spacing"] = h
             mu_a_set.attrs["spacing_units"] = "mm"
             mu_a_set.attrs["fov"] = fov
             mu_a_set.attrs["imaging_time"] = np.array(imaging_times)
             mu_a_set.attrs["wavelength"] = np.array(wavelengths)

            


            
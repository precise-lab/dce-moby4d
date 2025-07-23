import dolfin as dl
import ufl
import numpy as np

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append("../")
import moby


sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp
import scipy.io as io


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Compute SO2', fromfile_prefix_chars='@')
    parser.add_argument('-f', '--fname', default = "/workspace/shared_data/Moby_multi_wave/mesh/")
    parser.add_argument('-o', '--output', default = "/workspace/shared_data/Moby_multi_wave/masks/")
    parser.add_argument('--pzmin', default = 93, type = int)
    parser.add_argument('--pzmax', default = 107, type = int)
    parser.add_argument('--start_frame',
                        default=0,
                        type=int,
                        help = "First frame index")
    parser.add_argument('--end_frame',
                        default=200,
                        type=int,
                        help = "End frame index")

    cl_args = parser.parse_args()

    comm = dl.MPI.comm_world
    tissueComposition = moby.TissueComposition.create()
    labels = tissueComposition.tissue2label


    comm = dl.MPI.comm_world
    rank  = comm.rank


    for i in range(cl_args.start_frame, cl_args.end_frame):
        print(i)
        mesh = dl.Mesh(dl.MPI.comm_world)
        if i > 0:
            with dl.XDMFFile(cl_args.fname + f"moby_mesh{i}.xdmf") as fid:
                fid.read(mesh)
                geo_dim = mesh.geometry().dim()
                c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
                fid.read(c_labels, "c_labels")
        else:
            with dl.XDMFFile(cl_args.fname + "moby_mesh_structured.xdmf") as fid:
                fid.read(mesh)
                geo_dim = mesh.geometry().dim()
                c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
                fid.read(c_labels, "c_labels")
        dx = dl.Measure("dx", subdomain_data=c_labels)

        #Vh = dl.FunctionSpace(mesh, 'CG', 1)  
        Vh = dl.FunctionSpace(mesh, "DG", 0)
        m_trial, m_test = dl.TrialFunction(Vh), dl.TestFunction(Vh)
        varf_m = m_trial*m_test*dx

        rhs =  dl.Constant(0)*m_test*dx + dl.Constant(1)*m_test*dx(labels['artery']) + dl.Constant(2)*m_test*dx(labels['tumor']) + dl.Constant(3)*m_test*dx(labels['tumor_core'])
        A, b = dl.assemble_system(varf_m, rhs, [])
        Mask = dl.Function(Vh, name = 'mu_sp')
        dl.solve(A, Mask.vector(), b, 'cg', 'jacobi')

        
        h = 0.15
        z_fov = [cl_args.pzmin*h, cl_args.pzmax*h]
        Nx = 248
        xi = np.linspace(.5, (Nx-.5), Nx)*h
        yi = np.linspace(.5, (Nx-.5), Nx)*h
        zi = np.arange(z_fov[0]+.5*h, z_fov[1], h)
        XX, YY, ZZ = np.meshgrid(xi, yi, zi)
        points = np.hstack([np.reshape(xyz, (xyz.size,1)) for xyz in [XX,YY,ZZ] ])

        B = hp.assemblePointwiseObservation(Vh, points)
        mask_np = (B*Mask.vector()).gather_on_zero()

        if rank == 0:
            io.savemat(cl_args.output + f"mask_{i}.mat", {'mask': np.reshape(mask_np, XX.shape).astype(int)}, do_compression=True)



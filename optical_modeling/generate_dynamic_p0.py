'''
Created on May 10, 2022

@author: uvilla
'''


import json
import argparse
import configparser

import dolfin as dl
import ufl
import mpi4py
import numpy as np

import scipy.io as sio

from timeit import default_timer as timer
import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../../hippylib-public") )
import hippylib as hp


def convert_to_cartesian_parallel(comm, p0, N, args):

    rank  = comm.rank

    Nx, Ny, Nz = N
    xi = np.linspace(.5, (Nx-.5), Nx)*args['voxel_size_x']
    yi = np.linspace(.5, (Ny-.5), Ny)*args['voxel_size_y']
    zi = np.linspace(.5, (Nz-.5), Nz)*args['voxel_size_z']

    XX, YY, ZZ = np.meshgrid(xi, yi, zi)

    points = np.hstack([np.reshape(xyz, (xyz.size,1)) for xyz in [XX,YY,ZZ] ])

    Vh = p0.function_space()
    B = hp.assemblePointwiseObservation(Vh, points)

    p0_np = (B*p0.vector()).gather_on_zero()

    if rank == 0:
        sio.savemat(args['p0_name']+'.mat', {'p0': np.reshape(p0_np, XX.shape).astype(np.float32)}, do_compression=True)




def compute_frame(comm, args):
    
    hx, hy, hz = args['voxel_size_x'], args['voxel_size_y'],args['voxel_size_z']
    
    geo_dim = 3
    
    mesh = dl.Mesh(comm)
    with dl.XDMFFile(args['mesh_name']+".xdmf") as fid:        
        fid.read(mesh)
        c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
        c_labels.rename("c_labels", "c_labels")
        fid.read(c_labels, "c_labels")
    
    data_properties = sio.loadmat(args['properties_name']+".mat")
    mu_a_np = data_properties['mu_a']
    mu_a_exp = hp.NumpyScalarExpression3D()
    mu_a_exp.setData(mu_a_np, hx, hy, hz)
    
    mu_sp_np = data_properties['mu_sp']
    
    mu_tot = np.maximum( 3.*(mu_a_np+mu_sp_np), 1e-6)
    D_np = 1./ mu_tot
    D_exp = hp.NumpyScalarExpression3D()
    D_exp.setData(D_np, hx, hy, hz)
    
    z_min = 50 -.5*34
    z_max = 50 +.5*34
    illumination = dl.Expression(".25*(1.+std::tanh( x[2] - ZMIN))*(1.+std::tanh(ZMAX  - x[2]))", degree=1, ZMIN=z_min, ZMAX=z_max)
    
    Vh = dl.FunctionSpace(mesh, "CG", 1)
    DG0 = dl.FunctionSpace(mesh, "DG", 0)
    
    phi_trial = dl.TrialFunction(Vh)
    phi_test  = dl.TestFunction(Vh)
    
    mu_a = dl.interpolate(mu_a_exp, DG0)
    mu_a.rename("mu_a", "mu_a")
    D = dl.interpolate(D_exp, DG0)
    
    dx_diff   = ufl.Measure("dx", subdomain_data=c_labels, metadata={"quadrature_degree": 0})
    dx_lumped = ufl.dx(metadata={"quadrature_degree": 1, "representation":"quadrature"}, scheme='vertex')
    ds_lumped = ufl.ds(metadata={"quadrature_degree": 1, "representation":"quadrature"}, scheme='vertex')

    z_dir = dl.Constant((0.,0., 1.))
    Aform = D*ufl.inner(ufl.grad(phi_trial), ufl.grad(phi_test))*dx_diff   \
            -dl.Constant(0.99)*D*ufl.inner(ufl.dot(z_dir, ufl.grad(phi_trial)), ufl.dot(z_dir, ufl.grad(phi_test)))*dx_diff(1) \
            + ufl.inner(mu_a*phi_trial, phi_test)*dx_lumped \
            + ufl.inner(dl.Constant(0.5)*phi_trial, phi_test)*ds_lumped
            
    bform = ufl.inner(dl.Constant(0.5)*illumination, phi_test)*ds_lumped
    
    fluence = dl.Function(Vh, name="Fluence")
    
    A,b = dl.assemble_system(Aform, bform)
    Asolver = hp.PETScKrylovSolver(comm, "cg", "hypre_amg")
    Asolver.set_operator(A)
    
    Asolver.solve(fluence.vector(), b)

    DG1 = dl.FunctionSpace(mesh, "DG", 1)
    p0 = dl.project(mu_a*fluence, solver_type='cg', preconditioner_type='jacobi')
    p0.rename("p0", "p0")
        
    with dl.XDMFFile(comm, args['fluence_name']+".xdmf") as fid:
        fid.parameters["functions_share_mesh"] = True
        fid.parameters["rewrite_function_mesh"] = False
        fid.write(fluence,0)
        fid.write(mu_a, 0)
        fid.write(p0, 0)
        labels_np = data_properties['label_map'].astype(np.float64)
        labels_exp =  hp.NumpyScalarExpression3D()
        labels_exp.setData(labels_np, hx, hy, hz)
        labels = dl.interpolate(labels_exp, DG0)
        labels.rename("labels", "labels")
        fid.write(labels,0)


    convert_to_cartesian_parallel(comm, p0, mu_a_np.shape, args)
    
    dx = dl.Measure("dx", subdomain_data=c_labels)
    lesion_volume = dl.assemble( (dl.Constant(1.)+ dl.Constant(0.)*fluence)*dx(args['lesion_label']))
    lesion_signal = dl.assemble( mu_a*fluence*dx(args['lesion_label']))/lesion_volume
    
    liver_volume = dl.assemble( (dl.Constant(1.)+ dl.Constant(0.)*fluence)*dx(args['spleen_label']))
    liver_signal = dl.assemble( mu_a*fluence*dx(args['spleen_label']))/liver_volume
    
    return lesion_signal, liver_signal, lesion_volume, liver_volume
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate unstructured mesh', fromfile_prefix_chars='@')
    parser.add_argument('--config',
                        default='/workspace/shared_data/DCE-MOBY4D/vit_z1_r3/config.ini',
                        type=str,
                        help = "Configuration filename")
    parser.add_argument('--start_frame',
                        default=0,
                        type=int,
                        help = "First frame index")
    parser.add_argument('--end_frame',
                        default=-1,
                        type=int,
                        help = "End frame index")

    cl_args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(cl_args.config)

    args = {}
    args['properties']=os.path.join(config.get('path', 'root_folder'), 'phantom', config.get('path','fname_template'))
    args['mesh']=os.path.join(config.get('path', 'root_folder'), 'mesh', config.get('path','fname_template'))
    args['fluence']=os.path.join(config.get('path','root_folder'), 'fluence', config.get('path','fname_template'))
    args['p0']=os.path.join(config.get('path','root_folder'), 'p0', config.get('path','fname_template'))
    args['label_map'] = 'label_map'
    args['voxel_size_x'] = config.getfloat('grid','voxel_size_x')
    args['voxel_size_y'] = config.getfloat('grid','voxel_size_y')
    args['voxel_size_z'] = config.getfloat('grid','voxel_size_z')
    args['n_anatomy_frames'] = config.getint('grid','anatomy_frame_n')
    args['n_frames'] = config.getint('grid','frame_n')

    args['lesion_label'] = config.getint('labels','lesion')+1
    args['spleen_label'] = config.getint('labels','spleen')+1

    nstart = cl_args.start_frame
    nend   = cl_args.end_frame if cl_args.end_frame > 0 else args['n_frames']

    comm = dl.MPI.comm_world  #2018.1
    rank  = comm.rank
    nproc = comm.size

    dl.set_log_active(False)

    if rank == 0:
        f = open('average_signal_{0}_{1}.txt'.format(nstart, nend), 'w')
        
    if rank == 0:
        print("index lesion spleen time\n") 
    for index in range(nstart, nend, 1):
        args['mesh_name'] = args['mesh'].format( index% args['n_anatomy_frames'] +1)
        args['properties_name'] = args['properties'].format(index+1)
        args['fluence_name'] = args['fluence'].format(index+1)
        args['p0_name'] = args['p0'].format(index+1)
        tic = timer()
        lesion_signal, liver_signal, lesion_volume, liver_volume =  compute_frame(comm, args)
        toc = timer()
        if rank == 0:
            f.write("{0:4d} {1:1.5e} {2:1.5e} {3:1.5e} {4:1.5e}\n".format(
                index, lesion_signal, liver_signal, lesion_volume, liver_volume))
            f.flush()
            print("{0:4d} {1:1.5e} {2:1.5e} {3:1.5e}".format(index, lesion_signal, liver_signal, toc-tic))

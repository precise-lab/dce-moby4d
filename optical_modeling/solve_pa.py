import numpy as np
import dolfin as dl
import ufl
import scipy.io as io
import time

import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib") )
from hippylib import *

sys.path.append('../../qpact2/')
from qpact2 import *

import argparse




def solve_pa(lamb, inds):
    nx = 208
    s = 14.
    pixel_size = 2*s/nx

    theta_max = 12.5*np.pi/180.
    intensity = 2e4
    mua = 1e-3
    source_locations = 1e3*io.loadmat('source_locations.mat')['xy']
    source_locations = source_locations[:,[0,2]] 
    source_locations = source_locations@np.array([[0, 1], [-1, 0]])
    source_locations = source_locations + s
    source_locations = source_locations[[0,-1],:]
    nsources = source_locations.shape[0]

    abspath = '../../qpact2/applications/source_model/'
    source_directory = abspath  
    fname = "point_source.cpp"
    
    with open(os.path.join(source_directory, fname), "r") as cpp_file:
        cpp_code  = cpp_file.read()
        include_dirs = [".", source_directory]
    cpp_module = dl.compile_cpp_code(cpp_code, include_dirs=include_dirs)


    obs_grid_marks = np.linspace(0,2*s, nx+1)[:-1] + 1e-6
    x0, x1 = np.meshgrid(obs_grid_marks, obs_grid_marks)

    mesh = dl.RectangleMesh(dl.Point(0, 0), dl.Point(2*s, 2*s), nx, nx)
    mu_a = NumpyScalarExpression2D()
    mu_sp = NumpyScalarExpression2D()


    
    out_folder = '/workspace/shared_data/p0_multi_wave/p0_'+str(lamb) + '/'

    try:
        os.mkdir(out_folder)
        print("Making folder: "+ out_folder)
    except:
        print("Folder " + out_folder+ " already exists")


    for j in inds:
        print("     " + str(j))
        
        M = io.loadmat('/workspace/shared_data/Moby_multi_wave/phantom_'+ str(lamb) +'/moby' + '0'*(6-len(str(j+1))) + str(j+1) + '.mat')
        
        
        mask = dl.MeshFunction('size_t', mesh, 2)
        np_mask = (np.transpose(M['label_map'][20:-20,20:-20]) >0)
        numpy2MeshFunction(mesh, pixel_size*np.ones(2, dtype=np.float64), np_mask.astype(int), mask)
        sub_mesh = dl.SubMesh(mesh, mask, 1)


        Vh_phi = dl.FunctionSpace(sub_mesh, "CG", 1) 
        Vh_pa  = dl.FunctionSpace(sub_mesh, "CG", 1)
        Vh_m   = dl.VectorFunctionSpace(sub_mesh, "CG", 1)
        Vh = [[Vh_phi, Vh_pa], Vh_m, [Vh_phi, Vh_pa]]


        illumination = dl.CompiledExpression( cpp_module.BoundaryQ0Source(sub_mesh), degree=1 )
        
        for i in range(nsources):
            S = cpp_module.ConeBeam(2)
            S.set(source_locations[i,:].flatten(), (s - source_locations[i,:].flatten())/40, theta_max, intensity, mua)
            illumination.append(S)
            #print(( (s - source_locations[i,0])**2 + (s - source_locations[i,1])**2)**0.5 /40)

        transform = IdentityTransform()
        phi_handler = SP1(illumination, dl.ds, transform)
        pa_handler  = PAmodel(phi_handler, dl.Constant(0.), dl.Constant(1.))

        pde = SFSI_QPACT_Problem(Vh, phi_handler, pa_handler, input_safe=True)
        u = pde.generate_state()
        
        

        mu_a.setData(np.transpose(M['mu_a'][20:-20,20:-20]), pixel_size,pixel_size)
        Mu_a = dl.interpolate(mu_a, Vh_pa)
        mu_sp.setData(np.transpose(M['mu_sp'][20:-20,20:-20]),pixel_size,pixel_size)
        Mu_sp = dl.interpolate(mu_sp, Vh_pa)

        m = ufl.as_vector((Mu_a, Mu_sp))
        x = [u, dl.project(m,Vh_m).vector(), None]
        pde.solveFwd(u, x)

        obs_points = np.stack((x1[np_mask], x0[np_mask]), 1)

        B = assemblePointwiseObservation(Vh_pa, obs_points)
        out = dl.Vector(B.mpi_comm())

        B.init_vector(out, 0)
        B.mult(u.pa, out)

        p0 = np.zeros((nx,nx))
        p0[np_mask] = out.get_local()
        
        """plt.imshow(p0**0.5, vmin = 0, vmax = 0.2)
        plt.colorbar()
        plt.savefig('../system_example.png')
        assert False"""

        io.savemat(out_folder + 'p0_'+str(j+1)+'.mat',{'p0':p0})

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--lamb', help='wavelength to solve pa', required=True, type = int)
    parser.add_argument('-s','--start_ind', help='frame to start solving', required=True, type = int)
    parser.add_argument('-e','--end_ind', help='frame to end solving', required=True, type = int)
    args = parser.parse_args()
    print(args)
    

    solve_pa(args.lamb, range(args.start_ind,args.end_ind))




    
    
    







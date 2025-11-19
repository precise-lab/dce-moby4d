'''
Created on May 10, 2022

@author: uvilla
'''
'''
Created on May 23, 2019

@author: uvilla
'''
import dolfin as dl
import numpy as np

import argparse


from timeit import default_timer as timer

import sys
import os
sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib") )
sys.path.append("../")
import moby

import h5py

import hippylib as hp

import os
    
def generate_structured_mesh(input_fname, output_fname, z_bnds):    
    with h5py.File(input_fname, "r") as fid:
        data_set = fid["/phantom"]
        h = data_set.attrs["spacing"]
        data = data_set[:,:,:].astype(np.int32)

    geo_dim = len(data.shape)
    assert( geo_dim == 3) 
    print("(nx, ny, nz) = ({0}, {1}, {2})".format(*data.shape) )

    z_bnds_mm = [z_bnds[0]*h[2], z_bnds[1]*h[2]]
    points = [dl.Point(0.0, 0.0, z_bnds_mm[0]),
              dl.Point(data.shape[0]*h[0], data.shape[1]*h[1],z_bnds_mm[1]) ]
    outer_mesh = dl.BoxMesh.create(points, data.shape, dl.cpp.mesh.CellType.Type.tetrahedron)
    mask = dl.MeshFunction('size_t', outer_mesh, 3)
    np_mask = data > 0
    hp.numpy2MeshFunction(outer_mesh, h, np_mask.astype(int), mask)
    mesh = dl.SubMesh(outer_mesh, mask, 1)
    c_labels = dl.MeshFunction('size_t', mesh, 3)
    c_labels.rename("c_labels", "c_labels")
    hp.numpy2MeshFunction(mesh, h, data, c_labels)
    
    print('Number of elements: ', mesh.num_cells())
    print('Number of vertices: ', mesh.num_vertices())
        
    with dl.XDMFFile(output_fname) as fid:
        fid.write(mesh)
        fid.write(c_labels)

        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate unstructured mesh', fromfile_prefix_chars='@')
    parser.add_argument('-f', '--fname', default = "phantom_anatomy_0.h5")
    parser.add_argument('-o', '--output', default = "moby_mesh.xdmf")
    parser.add_argument('--zmin', default = 150, type = int)
    parser.add_argument('--zmax', default = 400, type = int)
    #parser.add_argument('-o', '--output', default = "/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh_structured.xdmf")

    cl_args = parser.parse_args()
    input_fname = cl_args.fname
    output_fname = cl_args.output
    z_bnds = [cl_args.zmin, cl_args.zmax]
    generate_structured_mesh(input_fname, output_fname, z_bnds)

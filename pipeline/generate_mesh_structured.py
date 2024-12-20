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
import pygalmesh
import meshio

from timeit import default_timer as timer

import sys
import os

sys.path.append("../")
import moby


sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib-public") )
import hippylib as hp

import os
    
def generate_structured_mesh(input_fname, output_fname, h):    
    data = np.load(input_fname)["phantom"]
    geo_dim = len(data.shape)
    assert( geo_dim == 3) 
    print("(nx, ny, nz) = ({0}, {1}, {2})".format(*data.shape) )

    z_data = np.max(data, axis=(0,1))

    z_indexes = np.where(z_data == moby.tissue2label['tumor'])[0]
    print(z_indexes[0], z_indexes[-1])
    z_index_min = z_indexes[0] - 100
    z_index_max = z_indexes[-1] + 100
    print("extended indexes", z_index_min, " ", z_index_max)

    data[:,:, 0:z_index_min] = moby.tissue2label["background"]
    data[:,:, z_index_max:] = moby.tissue2label["background"]

    outer_mesh = dl.BoxMesh(dl.Point(0.0, 0.0, 0.0), dl.Point(data.shape[0]*h[0], data.shape[1]*h[1],data.shape[2]*h[2]), *data.shape)
    mask = dl.MeshFunction('size_t', outer_mesh, 3)
    np_mask = data > 0
    hp.numpy2MeshFunction(outer_mesh, np.array(h), np_mask.astype(int), mask)
    mesh = dl.SubMesh(outer_mesh, mask, 1)
    c_labels = dl.MeshFunction('size_t', mesh, 3)
    c_labels.rename("c_labels", "c_labels")
    hp.numpy2MeshFunction(mesh, np.array(h), data, c_labels)
    
    print('Number of elements: ', mesh.num_cells())
    print('Number of vertices: ', mesh.num_vertices())
    
    
    with dl.XDMFFile(output_fname) as fid:
        fid.write(mesh)
        fid.write(c_labels)
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate unstructured mesh', fromfile_prefix_chars='@')
    parser.add_argument('-hx', default = 0.15, help="Grid size x-direction")
    parser.add_argument('-hy', default =0.15, help="Grid size y-direction")
    parser.add_argument('-hz', default = 0.15, help="Grid size z-direction")
    parser.add_argument('-f', '--fname', default = "../dce-moby-lfs/moby_phantom.npz")
    parser.add_argument('-o', '--output', default = "moby_mesh.xdmf")
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help = "Activate verbose output of pygalmesh")


    cl_args = parser.parse_args()
    input_fname = cl_args.fname
    output_fname = cl_args.output
    h = [cl_args.hx, cl_args.hy, cl_args.hz]
    verbose = cl_args.verbose
    generate_structured_mesh(input_fname, output_fname, h)

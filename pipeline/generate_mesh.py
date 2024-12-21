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
sys.path.append("../")
import moby

import os



def generate_mesh(data, h, cell_sizes_map, verbose):
    
    dd = data.astype(dtype=np.uint16)

    mesh = pygalmesh.generate_from_array(dd, h,lloyd=False, odt=False, max_cell_circumradius=cell_sizes_map,
                                         max_facet_distance=1.25*h[0], verbose=verbose)

    
    dd_unique = np.unique(dd[:]).astype(np.uint32)
    cells = mesh.get_cells_type('tetra')
    old_labels = mesh.get_cell_data("medit:ref", 'tetra')
    labels = dd_unique[old_labels]
   
    print(labels.shape, cells.shape)

    mesh.cells = [meshio.CellBlock('tetra', cells)] 
    mesh.cell_data["c_labels"] = [labels]

    fname='tmp{0}'.format(np.random.randint(10000,100000-1))
    meshio.xdmf.write(fname+".xdmf", mesh)
    
    dl_mesh = dl.Mesh()
    with dl.XDMFFile(fname+".xdmf") as fid:
        fid.read(dl_mesh)
        geo_dim = dl_mesh.geometry().dim()
        c_labels = dl.MeshFunction('size_t', dl_mesh, geo_dim)
        c_labels.rename("c_labels", "c_labels")
        fid.read(c_labels, "c_labels")
        
    os.remove(fname+".xdmf")
    os.remove(fname+".h5")
    

    return dl_mesh, c_labels

    
def generate_unstructured_mesh(input_fname, output_fname, h, verbose):    
    data = np.load(input_fname)["phantom"]
    geo_dim = len(data.shape)
    assert( geo_dim == 3) 
    print("(nx, ny, nz) = ({0}, {1}, {2})".format(*data.shape) )

    tissueComposition = moby.TissueComposition.create()
    tissue2label = tissueComposition.tissue2label
    data = tissueComposition.gLabelMap(data)

    z_data = np.max(data, axis=(0,1))

    z_indexes = np.where(z_data == tissue2label['tumor'])[0]
    print(z_indexes[0], z_indexes[-1])
    z_index_min = z_indexes[0] - 100
    z_index_max = z_indexes[-1] + 100
    print("extended indexes", z_index_min, " ", z_index_max)

    data[:,:, 0:z_index_min] = tissue2label["background"]
    data[:,:, z_index_max:] = tissue2label["background"]

    cell_sizes_map = {}
    cell_sizes_map['default'] = .5 #mm
    cell_sizes_map[tissue2label["artery"]] = 0.15 #mm
    cell_sizes_map[tissue2label["vein"]] = 0.15 #mm
    cell_sizes_map[tissue2label["tumor"]] = 0.15
    cell_sizes_map[tissue2label["liver"]] = 0.25
    cell_sizes_map[tissue2label["gall_bladder"]] = 0.25

    start = timer()
    mesh, c_labels = generate_mesh(data, h, cell_sizes_map, verbose)
    end = timer()
    print('Generate Mesh: Elapsed time: {}'.format(end-start))
    
    print('Number of elements: ',  mesh.num_cells())
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
    generate_unstructured_mesh(input_fname, output_fname, h, verbose)

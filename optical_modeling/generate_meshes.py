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
import configparser
import scipy.io as sio
import pygalmesh
import meshio

from timeit import default_timer as timer

import sys
import os



def generate_mesh(data, h):
    

    dd = data.astype(dtype=np.uint16)

    cell_sizes_map = {}
    cell_sizes_map['default'] = 0.5
        
    

    mesh = pygalmesh.generate_from_array(dd, h, max_cell_circumradius=cell_sizes_map,
                                         max_facet_distance=.5*h[0], verbose=True)

    
    dd_unique = np.unique(dd[:]).astype(np.uint32)
    #try:
    #    mesh.remove_lower_dimensional_cells()
    #except:
    #    mesh.prune()
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

    
def generate_unstructured_mesh(args, index):    
    phantom = sio.loadmat(args['name'])

    data = phantom[args['label_map']]
    geo_dim = len(data.shape)
    h = [args['voxel_size_x'], args['voxel_size_y'], args['voxel_size_z']]
    
    Nx = data.shape[0]
    Ny = data.shape[1]
    Nz = data.shape[2]
    nxi = np.arange(0, Nx)
    nyi = np.arange(0, Ny)
    nzi = np.arange(0, Nz)
    
    print("(nx, ny) = ({0}, {1})".format(data.shape[0], data.shape[1]) )
    
    dxi = 2.*(nxi-.5*Nx)/Nx
    dyi = 2.*(nyi-.5*Ny)/Ny
    dzi = 2.*(nzi-.5*Nz)/Nz

    dxi, dyi, dzi = np.meshgrid(dxi, dyi, dzi)
    
    data[dxi*dxi+dyi*dyi >= 1.] = -1
    data = data+1
     
    start = timer()
    mesh, c_labels = generate_mesh(data, h)
    end = timer()
    print('Generate Mesh: Elapsed time: {}'.format(end-start))
    
    print('Number of elements: ',  mesh.num_cells())
    print('Number of vertices: ', mesh.num_vertices())
    
    
    with dl.XDMFFile(args['out_name']+".xdmf") as fid:
        fid.write(mesh)
        fid.write(c_labels)
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate unstructured mesh', fromfile_prefix_chars='@')
    parser.add_argument('--config',
                        default='/workspace/shared_data/DCE-MOBY4D/output1/config.ini',
                        type=str,
                        help = "Configuration filename")

    cl_args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(cl_args.config)

    #root_folder = /workspace/shared_data/DCE-MOBY4D/output1
    #fname_template= moby{0:06d}
    #frame_n = 13680
    #anatomy_frame_n = 287
    #voxel_size_x = 0.150000
    #voxel_size_y = 0.150000
    #voxel_size_z = 0.150000
    #frame_rate = 0.050000

    args = {}
    args['input']=os.path.join(config.get('path', 'root_folder'), config.get('path','fname_template'))
    args['output']=os.path.join(config.get('path','root_folder'), '../mesh', config.get('path','fname_template'))
    args['label_map'] = 'label_map'
    args['voxel_size_x'] = config.getfloat('grid','voxel_size_x')
    args['voxel_size_y'] = config.getfloat('grid','voxel_size_y')
    args['voxel_size_z'] = config.getfloat('grid','voxel_size_z')
    args['n_frames'] = config.getint('grid','anatomy_frame_n')
    
    
    for index in range(args['n_frames']):
        args['name'] = args['input'].format(index+1)
        args['out_name'] = args['output'].format(index+1)
        print(args['name'])
        generate_unstructured_mesh(args, index)

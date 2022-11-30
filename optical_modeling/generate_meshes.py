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

import json
import argparse
import configparser
import scipy.io as sio
import pygalmesh
import meshio

from timeit import default_timer as timer

import sys
import os



def generate_mesh(data, args):
    
    h = [args['voxel_size_x'], args['voxel_size_y'], args['voxel_size_z']]
    dd = data.astype(dtype=np.uint16)

    cell_sizes_map = {}
    cell_sizes_map['default'] = .51 #mm
    cell_sizes_map[args['spleen_label']+1] = 0.25 #mm
    cell_sizes_map[args['lesion_label']+1] = 0.25 #mm
        
    

    mesh = pygalmesh.generate_from_array(dd, h, max_cell_circumradius=cell_sizes_map,
                                         max_facet_distance=1.1*h[0], verbose=args['verbose'])

    
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

    
def generate_unstructured_mesh(args, index):    
    phantom = sio.loadmat(args['name'])

    data = phantom[args['label_map']]
    geo_dim = len(data.shape)

    default_label = args['intestin_label'][0]
    for l in args['intestin_label']:
        data[data==l] = default_label

    default_label = args['brain_label'][0]
    for l in range(default_label, args['brain_label'][1]+1):
        data[data==l] = default_label
    
    
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
    mesh, c_labels = generate_mesh(data, args)
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
                        default='/workspace/shared_data/DCE-MOBY4D/vit_z1_r3/config.ini',
                        type=str,
                        help = "Configuration filename")
    parser.add_argument('--start_frame',
                        default=1,
                        type=int,
                        help = "First frame to be meshed")
    parser.add_argument('--last_frame',
                        default=-1,
                        type=int,
                        help = "Last frame to be meshed")
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help = "Activate verbose output of pygalmesh")
    parser.add_argument('-s', '--skip',
                        action='store_false',
                        help = "Skip if file exists")

    cl_args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(cl_args.config)

    args = {}
    args['input']=os.path.join(config.get('path', 'root_folder'), 'phantom', config.get('path','fname_template'))
    args['output']=os.path.join(config.get('path','root_folder'), 'mesh', config.get('path','fname_template'))
    args['label_map'] = 'label_map'
    args['voxel_size_x'] = config.getfloat('grid','voxel_size_x')
    args['voxel_size_y'] = config.getfloat('grid','voxel_size_y')
    args['voxel_size_z'] = config.getfloat('grid','voxel_size_z')
    args['n_frames'] = config.getint('grid','anatomy_frame_n')
    args['spleen_label'] = config.getint('labels','spleen')
    args['lesion_label'] = config.getint('labels','lesion')
    args['intestin_label'] = json.loads(config.get('labels','intestin'))
    args['brain_label'] = json.loads(config.get('labels','brain'))
    
    last_frame = cl_args.last_frame if (cl_args.last_frame > 0) else args['n_frames'] 
    for index in range(cl_args.start_frame, last_frame+1):
        args['name'] = args['input'].format(index)
        args['out_name'] = args['output'].format(index)
        args['verbose'] = cl_args.verbose
        print(args['name'])
        if os.path.isfile(args['out_name']+'.xdmf') and cl_args.skip:
            print("Skipping", args['out_name'])
        else:
            generate_unstructured_mesh(args, index)

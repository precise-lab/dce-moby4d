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
    
def move_mesh(ref_mesh, ref_labels, ref_Vh, disp_array, h):
    mesh = dl.Mesh(ref_mesh)

    for i in range(3):
        disp_array[:,:,:, i] *= 3.*h[i]

    disp_exp = hp.NumpyVectorExpression3D(3)
    disp_exp.setData(disp_array, *h)
    d = dl.interpolate(disp_exp, ref_Vh)
    d_np = np.abs( d.vector().gather_on_zero() )
    d_np = np.reshape(d_np, (int(d_np.shape[0]/3), 3))
    nd_np = np.sqrt(np.sum(d_np*d_np, axis=1))
    if nd_np.size > 0:
        print("-->",  np.max(nd_np), " ", np.max(d_np[:,0]), " ", np.max(d_np[:,1]), " ", np.max(d_np[:,2]))

    dl.ALE.move(mesh, d)

    c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
    c_labels.array()[:] = ref_labels.array()[:]
    c_labels.rename("c_labels", "c_labels")

    DG0 = dl.FunctionSpace(mesh, "DG", 0)
    test = dl.TestFunction(DG0)
    vol = dl.assemble(test*ufl.dx)
    vol_np = vol.gather_on_zero()
    if vol_np.size > 0:
        print( np.min(vol_np), " ", np.max(vol_np) )

    return mesh, c_labels, d

if __name__ == "__main__":
    nframes = 200
    h = [0.15, 0.15, 0.15]
    ref_mesh = dl.Mesh()
    #with dl.XDMFFile("moby_mesh.xdmf") as fid:
    #with dl.XDMFFile("/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh.xdmf") as fid:
    with dl.XDMFFile("/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh_structured.xdmf") as fid:
   
        fid.read(ref_mesh)
        geo_dim = ref_mesh.geometry().dim()
        ref_labels = dl.MeshFunction('size_t', ref_mesh, geo_dim)
        ref_labels.rename("c_labels", "c_labels")
        fid.read(ref_labels, "c_labels")
    print(ref_mesh.num_vertices())
    ref_Vh = dl.VectorFunctionSpace(ref_mesh, "CG", 1)

    fid_out = dl.XDMFFile("displacement.xdmf")
    fid_out.parameters["functions_share_mesh"] = True
    fid_out.parameters["rewrite_function_mesh"] = False
    zmin = 175
    zmax = 375
    #fid_out.write(ref_labels)

    for i in range(84,nframes):
        #disp_array = np.load(f"../dce-moby-lfs/motion/motion{i}.npz")["motion"]
        disp_array = np.load(f"/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_array/motion{i}.npz")["motion"]
        disp_array = disp_array[:,:,zmin:zmax, :]
        mesh, c_labels, d = move_mesh(ref_mesh, ref_labels, ref_Vh, disp_array, h)
        #with dl.XDMFFile(f"moving_mesh/moby_mesh{i}.xdmf") as fid:
        with dl.XDMFFile(f"/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh{i}.xdmf") as fid:
        
            fid.write(mesh)
            fid.write(c_labels)

        fid_out.write(d, i)
        print(i)


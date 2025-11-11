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
    parser.add_argument('-w', '--wavelength', default = 800, type = float)
    parser.add_argument('-m', '--mesh', default = "moby_mesh.xdmf")
    
    args = parser.parse_args()

    comm = dl.MPI.comm_world

    wavelength = args.wavelength
    imaging_time = 0.
    fov = np.array([[0., 37.2], [0., 37.2], [26.25, 56.25]])
    N = [248, 248, 200]

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
    Vh_p0 = dl.FunctionSpace(mesh, "DG", 1)

    resampler = moby.CartesianGridResampler(Vh_p0, fov, N, subsampling=1)

    b_box = get_boundingbox(Vh_phi)
    if 0==comm.rank:
        print(b_box)
        
    center = .5*(b_box[:,0] + b_box[:,1])
    
    slit_height = 20
    n_sources = 20
    sources_distance = 40
    sources_angles = np.linspace(0, 2.*np.pi, n_sources, endpoint=False)
    theta_max = 9.5*np.pi/180.
    intensity = 100.
    source_mua = chromophoresDec.get_pure_mu_a(wavelength, [moby.Chromophore.WATER])[moby.Chromophore.WATER]
    if 0==comm.rank:
         print("Water mu_a: ", source_mua)


    illumination = dl.CompiledExpression( moby.BoundaryQ0Source(mesh), degree=1 )
    for isource in range(n_sources):
        source_direction = np.array([np.cos(sources_angles[isource]), np.sin(sources_angles[isource]), 0.])
        source_location = center - sources_distance*source_direction
        S = moby.SlitBeam(3)
        S.set(source_location.flatten(), source_direction.flatten(), slit_height,theta_max, intensity, source_mua)
        illumination.append(S)
            
    dx_diff   = ufl.Measure("dx", subdomain_data=c_labels, metadata={"quadrature_degree": 0})
    dx_lumped = ufl.dx(metadata={"quadrature_degree": 1, "representation":"quadrature"}, scheme='vertex')
    ds_lumped = ufl.ds(metadata={"quadrature_degree": 1, "representation":"quadrature"}, scheme='vertex')

    femPhantom = moby.FEMPhantom(Vh_m, Vh_phi, c_labels, tissueComposition, chromophoresDec, pkModel)
    sat_map = {"artery": 0.98, "vein": 0.7, "tumor": .371, "tumor_core": 0}
    so2 = femPhantom.compute_oxygen_saturation(sat_map )

    mu_a = femPhantom.compute_mu_a(wavelength, imaging_time, so2)
    mu_sp = femPhantom.compute_mu_sp(wavelength)
    D = 1./(3.*(mu_a + mu_sp))

    phi_trial, phi_test = dl.TrialFunction(Vh_phi), dl.TestFunction(Vh_phi)
    Aform = D*ufl.inner(ufl.grad(phi_trial), ufl.grad(phi_test))*dx_diff  \
                    + ufl.inner(mu_a*phi_trial, phi_test)*dx_lumped \
                    + ufl.inner(dl.Constant(0.5)*phi_trial, phi_test)*ds_lumped
    bform = ufl.inner(dl.Constant(0.5)*illumination, phi_test)*ds_lumped

                    
    fluence = dl.Function(Vh_phi, name="Fluence")
    A,b = dl.assemble_system(Aform, bform)
    Asolver = hp.PETScKrylovSolver(comm, "cg", "hypre_amg")
    Asolver.set_operator(A)
    Asolver.solve(fluence.vector(), b)

    max_fluence = fluence.vector().norm("linf")
    if 0==comm.rank:
        print("Max fluence: ", max_fluence)


    p0 = dl.project(mu_a*fluence, Vh_p0, solver_type='cg', preconditioner_type='jacobi')
    p0.rename("p0", "p0")

    with dl.XDMFFile(comm, "out.xdmf") as fid:
            fid.parameters["functions_share_mesh"] = True
            fid.parameters["rewrite_function_mesh"] = False
            fid.write(fluence,0)
            fid.write(mu_a, 0)
            fid.write(mu_sp, 0)
            fid.write(p0, 0)

    p0_np = resampler(p0)
    if 0 == comm.rank:
        pl = pv.Plotter()
        pl.add_volume(p0_np, cmap='viridis', opacity="sigmoid")
        pl.show(screenshot='p0.png')

    if 0 == comm.rank:
        with h5py.File("p0.h5", "w") as fid:
             d_set = fid.create_dataset("p0", data=p0_np, compression="lzf")
             h = (fov[:,1] - fov[:,0])/np.array(N)
             d_set.attrs["spacing"] = h
             d_set.attrs["spacing_units"] = "mm"
             d_set.attrs["fov"] = fov
             d_set.attrs["time"] = imaging_time
             d_set.attrs["wavelength"] = wavelength
            


            
import dolfin as dl
try:
    import ufl_legacy as ufl
except:
    import ufl
import numpy as np
import h5py

import inspect
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
    parser.add_argument('-p', '--phase', default=0, type=int)
    parser.add_argument('-m', '--mesh', default = "moby_mesh.xdmf")
    parser.add_argument('-d', '--disp', default = 'motion_vectors.h5')
    
    args = parser.parse_args()

    comm = dl.MPI.comm_world

    wavelength = args.wavelength
    resp_phase = args.phase

    fov = np.array([[0., 37.2], [0., 37.2], [26.25, 56.25]])
    N = [248, 248, 200]
    start_time = -36.
    end_time   = 684.01

    tissueComposition = moby.TissueComposition.create()
    chromophores = [moby.Chromophore.HB, moby.Chromophore.HBO2, moby.Chromophore.WATER, moby.Chromophore.CA]

    time = np.arange(-36., 684.01, 0.05)
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

    dx_labels = ufl.Measure("dx", domain = mesh, subdomain_data=c_labels)

    mover = moby.MeshMover(mesh, args.disp, verbose=True)
    dt = mover.dt
    number_resp_phases = mover.nframes

    if resp_phase > 0:
        mover.move2(dt*resp_phase, dx_labels)

    Vh_phi = dl.FunctionSpace(mesh, "CG", 1)
    Vh_m  = dl.FunctionSpace(mesh, "DG", 0)
    Vh_p0 = dl.FunctionSpace(mesh, "DG", 1)

    resampler_p0 = moby.CartesianGridResampler(Vh_p0, fov, N, subsampling=1)
    resampler_ind = moby.CartesianGridResampler(Vh_m, fov, N, subsampling=1)

    b_box = get_boundingbox(Vh_phi)
    if 0==comm.rank:
        print("Bounding box: ", b_box)
        
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
    mu_sp = femPhantom.compute_mu_sp(wavelength)
 
    indicator_np = {}
    for key in ["artery", "tumor", "tumor_core", "liver"]:
        indicator = femPhantom.computeIndicatorFunction(key)
        indicator_np[key] = resampler_ind(indicator)

    if 0 == comm.rank:
        fname = f"indicator/indicator_{resp_phase:05d}.h5"
        spacing_indicator = (fov[:,1] - fov[:,0])/np.array(N)
        with h5py.File(fname, "w") as fid:
            for key in indicator_np:
                i_set = fid.create_dataset(key, data=indicator_np[key], compression="lzf")
                i_set.attrs["spacing"] = spacing_indicator
                i_set.attrs["spacing_units"] = "mm"
                i_set.attrs["fov"] = fov
                i_set.attrs["phase"] = resp_phase
                i_set.attrs["num_phases"] = number_resp_phases
                i_set.attrs["dt"] = dt

    imaging_times = np.arange(start_time+resp_phase*dt, end_time, dt*number_resp_phases)

    if 0 == comm.rank:
        fname = f"p0/p0_{resp_phase:05d}.h5"
        fid = h5py.File(fname, "w")
        spacing = (fov[:,1] - fov[:,0])/np.array(N)
        p0_size = tuple(N) + imaging_times.shape
        p0_chunck_size = tuple([ni//4 for ni in N]) + (1,)
        print("p0_size", p0_size)
        print("p0_chunck_size", p0_chunck_size)
        p0_set = fid.create_dataset("p0", p0_size, chunks=p0_chunck_size, compression="lzf", dtype=np.float32)
        p0_set.attrs["spacing"] = spacing
        p0_set.attrs["spacing_units"] = "mm"
        p0_set.attrs["fov"] = fov
        p0_set.attrs["resp_phase"] = resp_phase
        p0_set.attrs["wavelength"] = wavelength

        time_set = fid.create_dataset("time", data = imaging_times.astype(np.float32), compression="lzf")
        time_set.attrs["units"] = "s"


    # Preallocate solution vector, rhs, and lhs matrix for solving the DA
    fluence = dl.Function(Vh_phi, name="Fluence")
    A = dl.PETScMatrix()
    b = dl.PETScVector()
    Asolver = hp.PETScKrylovSolver(comm, "cg", "hypre_amg")
    Asolver.parameters["error_on_nonconvergence"] = False

    # Preallocate solution vector, rhs, and lhs matrix for computing p0
    p0 = dl.Function(Vh_p0, name="p0")
    p0_trial, p0_test = dl.TrialFunction(Vh_p0), dl.TestFunction(Vh_p0)
    m_po_varf = ufl.inner(p0_trial, p0_test) * ufl.dx
    Mp0 = dl.assemble(m_po_varf)
    Mp0_solver = hp.PETScKrylovSolver(comm, "cg", "jacobi")
    Mp0_solver.set_operator(Mp0)
    p0_rhs = dl.PETScVector()

    for index, imaging_time in enumerate(imaging_times[:]):
        tic = timer()
        mu_a = femPhantom.compute_mu_a(wavelength, imaging_time, so2)
        D = 1./(3.*(mu_a + mu_sp))
        
        phi_trial, phi_test = dl.TrialFunction(Vh_phi), dl.TestFunction(Vh_phi)
        Aform = D*ufl.inner(ufl.grad(phi_trial), ufl.grad(phi_test))*dx_diff  \
                    + ufl.inner(mu_a*phi_trial, phi_test)*dx_lumped \
                    + ufl.inner(dl.Constant(0.5)*phi_trial, phi_test)*ds_lumped
         
        bform = ufl.inner(dl.Constant(0.5)*illumination, phi_test)*ds_lumped

        dl.assemble_system(Aform, bform, A_tensor = A, b_tensor = b)
        Asolver.set_operator(A)
        Asolver.solve(fluence.vector(), b)
        reason = Asolver.ksp().getConvergedReason()
        
        if reason < 0:
            if comm.rank == 0:
                print(f"{__file__}:{inspect.currentframe().f_lineno}: Solver failed to converge! PETSc Reason code: {reason}")
                fid.close()
            sys.exit( reason )
            comm.Abort( reason  )

        dl.assemble( fluence*mu_a*p0_test*ufl.dx, tensor=p0_rhs)
        Mp0_solver.solve(p0.vector(), p0_rhs)
        
        
        p0_np = resampler_p0(p0)
        toc = timer()
        tictoc = toc - tic

        max_fluence = fluence.vector().norm("linf")
        if 0==comm.rank:
            print(f"Index: {index}; Time: {imaging_time}; Max fluence: {max_fluence}; Elapsed time: {tictoc:.6f}sec")

        if 0 == comm.rank:
            p0_set[:,:,:, index] = p0_np
            fid.flush()

            


            

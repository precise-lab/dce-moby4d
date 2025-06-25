import dolfin as dl
import ufl
import numpy as np
import scipy.io as io

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append("../")
import moby
import time


sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Compute SO2', fromfile_prefix_chars='@')
    parser.add_argument('-w', '--wavelength', default = 730, type = int)
    parser.add_argument('-p', '--mu_a', default = "/workspace/shared_data/Moby_multi_wave/mu_a_730/")
    parser.add_argument('-s', '--so2', default = "/workspace/shared_data/Moby_multi_wave/so2/")
    parser.add_argument('-m', '--mesh', default = "/workspace/shared_data/Moby_multi_wave/mesh/")
    parser.add_argument('--wb', default = 0, type = int)
    parser.add_argument('--nw', default = 8, type = int)

    parser.add_argument('--start_frame',
                        default=0,
                        type=int,
                        help = "First frame index")
    parser.add_argument('--end_frame',
                        default=9600,
                        type=int,
                        help = "End frame index")
    args = parser.parse_args()


    wavelength_batch = args.wb
    n_wave_lengths = args.nw

    comm = dl.MPI.comm_world
    rank  = comm.rank
    constants = io.loadmat('../properties2/constants.mat')
    constant_icg = io.loadmat('../properties2/constant_icg.mat')
    tissueComposition = moby.TissueComposition.create()
    labels = tissueComposition.tissue2label
    opt_prop = io.loadmat('../properties2/opt_prop.mat')

    A = 120 * 1e-6  
    Ap = 3
    B = 4.34/60
    C = 0.8*1e-6
    D = 1/60
    E = 0.07/60

    Ktrans  = 0.3/60
    kep     = 0.75/60    

    dt = 0.05
    startTime = -30
    endTime = 450
    t = np.arange( int((endTime - startTime)/dt))*dt + startTime
    t_pos = t[t> 0]
    ca_tp = np.zeros(t.shape)
    ca_tp[t>0] = A*((t_pos/60)**Ap)*np.exp(-B*t_pos) + C*(1 - np.exp(-D*t_pos))*np.exp(-E*t_pos)
    c_perf = np.zeros(t.shape)
    c_perf[t> 0] = (Ktrans*np.convolve(np.exp(-kep*t_pos), ca_tp[t>0])*dt)[:t_pos.shape[0]]
    del t_pos    

    mu_a_b_oxy = np.log(10)*tissueComposition.c_thb_b*constants['e_hbo2']
    mu_a_b_deoxy = np.log(10)*tissueComposition.c_thb_b*constants['e_hb']
    mu_a_w = constants['mu_a_w']

    g = opt_prop['g']
    mu_s_ref = opt_prop['mu_s_ref']
    wavelength_ref = opt_prop['wavelength_ref']
    wavelengths = constants['wavelength']


    f_ICG_in_blood = 50#;%50e-6/2e-3
    mu_a_ca = f_ICG_in_blood*np.log(10)*constant_icg['e_icg'][:,:1]

    nmesh = 199

    for mesh_it in range(nmesh):
        if rank == 0:
            print(f"Mesh {mesh_it}")
        mesh = dl.Mesh(comm)
        with dl.XDMFFile(args.mesh + f"moby_mesh{mesh_it%nmesh+1}.xdmf") as fid:
            fid.read(mesh)
            geo_dim = mesh.geometry().dim()
            c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
            fid.read(c_labels, "c_labels")
        dx = dl.Measure("dx", subdomain_data=c_labels, domain = mesh)
        Vols = [dl.assemble(dl.Constant(1.)*dx(labels[key])) for key in labels.keys()]

        Vols = [dl.assemble(dl.Constant(1.)*dx(labels[key])) for key in labels.keys()]
        Vh_m = dl.FunctionSpace(mesh, 'DG', 0)  
        m_trial, m_test = dl.TrialFunction(Vh_m), dl.TestFunction(Vh_m)

        varf_m = m_trial*m_test*dx
        rhs_musp =  dl.Constant(0.)*m_test*dx

        Vh_phi = dl.FunctionSpace(mesh, 'CG', 1)
        phi_trial, phi_test = dl.TrialFunction(Vh_phi), dl.TestFunction(Vh_phi)
        so2 = dl.Function(Vh_phi)
        with dl.XDMFFile(args.so2 + f"so2_{mesh_it%nmesh+1}.xdmf") as fid:
            fid.read_checkpoint(so2, "so2")

        h = 0.15
        
        z_fov = [93*h, 107*h]
        Nx = 248
        xi = np.linspace(.5, (Nx-.5), Nx)*h
        yi = np.linspace(.5, (Nx-.5), Nx)*h
        zi = np.arange(z_fov[0]+.5*h, z_fov[1], h)
        XX, YY, ZZ = np.meshgrid(xi, yi, zi)
        points = np.hstack([np.reshape(xyz, (xyz.size,1)) for xyz in [XX,YY,ZZ] ])

        B = hp.assemblePointwiseObservation(Vh_m, points)
            
        i = mesh_it #+ 16*nmesh
        while i < args.start_frame:
                i += nmesh
        while i < args.end_frame:
            if i%n_wave_lengths == wavelength_batch:
                #T1 = time.time()
                if rank == 0:
                    print(f"    Frame {i}")
                rhs_mua = dl.Constant(0.)*m_test*dx  

                for key, vol in zip(labels.keys(), Vols):
                    if vol > 0:
                        f_b = float(tissueComposition.volume_fractions[labels[key]][0])
                        f_w = float(tissueComposition.volume_fractions[labels[key]][1])
                        f_ca = f_b*ca_tp[i]
                        if key == 'tumor':
                            f_ca +=  (1. - f_b)*c_perf[i]
                        if key == 'tumor_core':
                            f_ca +=  0.5*(1. - f_b)*c_perf[i]
                        
                        s_coeff = mu_a_b_oxy[wavelengths == args.wavelength][0]*f_b\
                            - mu_a_b_deoxy[wavelengths == args.wavelength][0]*f_b
                        c_coeff = mu_a_b_deoxy[wavelengths == args.wavelength][0]*f_b \
                            + mu_a_w[wavelengths == args.wavelength][0]*f_w \
                            + mu_a_ca[wavelengths == args.wavelength][0]*f_ca
                        rhs_mua += (s_coeff*so2 + dl.Constant(c_coeff))*m_test*dx(labels[key])
                
                A, b = dl.assemble_system(varf_m, rhs_mua, [])
                Mu_a = dl.Function(Vh_m, name = 'mu_a')
                dl.solve(A, Mu_a.vector(), b, 'cg', 'jacobi')

                mu_a_np = (B*Mu_a.vector()).gather_on_zero()
                if rank == 0:
                    io.savemat(args.mu_a+ f"moby_{i}.mat", {'mu_a': np.reshape(mu_a_np, XX.shape).astype(np.float32),}, do_compression=True)
            i += nmesh
            


            
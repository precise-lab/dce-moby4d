import numpy as np
import dolfin as dl
import ufl
import scipy.io as io
import scipy.interpolate as interp
import scipy.sparse as sp
import time

from pylops import MatrixMult
from pylops.signalprocessing import Convolve1D, Convolve2D
from pylops.optimization.leastsquares import regularized_inversion

import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib") )
from hippylib import *

sys.path.append('../../qpact2/')
from qpact2 import *

import argparse

if __name__ == "__main__":
    
    imaging_op_folder = '/workspace/shared_data/msot_imaging_op/'
    p0_folder =  '/workspace/shared_data/Moby_multi_wave/'

    Nx = 248
    nx = 248
    
    ndet = 256
    ntm = 2030
    
    a = sp.load_npz("imaging_op.npz")
    imp_resp = np.load("imp_resp.npy")
    convolve_irf = Convolve1D((ndet, ntm), imp_resp, axis=1, offset=ntm // 2,
                                      dtype=np.float32)

    imaging_op = convolve_irf @ MatrixMult(a)

    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 1

    nT = 9600
    nt = nT//(ssf*nl)

    measurements = np.zeros((ndet, ntm, nl, nt))


    F = np.zeros((Nx,Nx,ssf))
    x = np.arange(nx)
    X = np.linspace(0,nx-1,Nx)
    for t in range(nt):
        print(t)
        for l in range(nl):
            for j in range(ssf):
                
                fid = l + j*ssf + t*ssf*nl

                fname = p0_folder + 'p0_' + str(lambs[l]) + '/moby_' + str(fid) + '.mat'
                f = io.loadmat(fname)['p0']
                #fi = interp.interp2d(x,x, f)
                #F[:,:,j] = fi(X,X)
                F[:,:,j] = np.mean(f,-1)

                #measurements[:,:,l,t] += imaging_op.matvec(F[:,:,j].flatten()).reshape((ndet,ntm))/ssf
                
                
           
            M = imaging_op.matmat(np.flip(F,0).reshape((Nx**2,ssf)))
            measurements[:,:,l,t] = np.mean(M,-1).reshape((ndet,ntm))

    np.save('acoustic_measurements.npy', measurements)
    io.savemat('acoustic_measurements.mat',{'measurements': measurements})






    



import numpy as np
import dolfin as dl
import ufl
import scipy.io as io
import scipy.interpolate as interp
import scipy.sparse as sp
import time

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
    p0_folder =  '/workspace/shared_data/p0_multi_wave/'

    Nx = 332
    nx = 208
    nm0 = 1249
    nm1 = 256
    vals = np.squeeze(io.loadmat(imaging_op_folder + 'imaging_vals.mat')['V'])
    ir = np.squeeze(io.loadmat(imaging_op_folder + 'imaging_rows.mat')['I']-1)
    jc = np.squeeze(io.loadmat(imaging_op_folder + 'imaging_cols.mat')['J']-1)
    a = sp.csc_matrix((vals, (ir, jc)), shape=(nm0*nm1, Nx**2))

    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 6

    nT = 9600
    nt = nT//(ssf*nl)

    measurements = np.zeros((nm0,nm1, nl, nt))


    F = np.zeros((Nx,Nx,ssf))
    x = np.arange(nx)
    X = np.linspace(0,nx-1,Nx)
    for t in range(nt):
        print(t)
        for l in range(nl):
            for j in range(ssf):
                
                fid = l + j*ssf + t*ssf*nl + 1

                fname = p0_folder + 'p0_' + str(lambs[l]) + '/p0_' + str(fid) + '.mat'
                f = io.loadmat(fname)['p0']
                fi = interp.interp2d(x,x, f)
                F[:,:,j] = fi(X,X)
           
            M = a@F.reshape((Nx**2,ssf),order = 'F')
            measurements[:,:,l,t] = np.mean(M,-1).reshape((nm0,nm1), order = 'F')
    io.savemat('acoustic_measurements.mat',{'measurements': measurements})






    



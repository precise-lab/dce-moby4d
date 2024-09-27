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
    
    p0_folder =  '/workspace/shared_data/Moby_multi_wave/'

    Nx = 332
    nx = 208


    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 6

    nT = 9600
    nt = nT//(ssf*nl)



    F = np.zeros((Nx,Nx, nt))
    G = np.zeros((Nx,Nx,nt))


    x = np.arange(nx)
    X = np.linspace(0,nx-1,Nx)
    for t in range(nt):
        print(t)
        for l in range(nl):
            for j in range(nl*ssf):
                fid = t*nl*ssf + j + 1
                fname = p0_folder + 'phantom_' + str(lambs[l]) + '/moby' + '0'*(6- len(str(fid))) +str(fid) + '.mat'
                m = io.loadmat(fname)['label_map'][20:-20, 20:-20]
                f = (m==9).astype(int)
                g = (m == 18).astype(int)
                fi = interp.interp2d(x,x, f)
                F[:,:,t] += fi(X,X)
                gi = interp.interp2d(x,x, g)
                G[:,:,t] += gi(X,X)
    F = (F > 0).astype(int)
    G = (G > 0).astype(int)
    io.savemat('subsample_mask.mat',{'mask': F})
    io.savemat('subsample_aif_mask.mat',{'aif_mask': G})
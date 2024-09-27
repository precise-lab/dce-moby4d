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
    
    p0_folder =  '/workspace/shared_data/p0_multi_wave/'

    Nx = 332
    nx = 208


    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 6

    nT = 9600
    nt = nT//(ssf*nl)



    F = np.zeros((Nx,Nx,nl, nt))


    x = np.arange(nx)
    X = np.linspace(0,nx-1,Nx)
    for t in range(nt):
        print(t)
        for l in range(nl):
            for j in range(nl*ssf):
                fid = t*nl*ssf + j + 1
                fname = p0_folder + 'p0_' + str(lambs[l]) + '/p0_' + str(fid) + '.mat'
                f = io.loadmat(fname)['p0']
                fi = interp.interp2d(x,x, f)
                F[:,:,l,t] += fi(X,X)/(nl*ssf)    
    io.savemat('subsample_p0.mat',{'p0': F})






    



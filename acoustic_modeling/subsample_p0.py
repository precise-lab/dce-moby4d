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

    Nx = 248
    nx = 248


    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 1

    nT = 9600
    nt = nT//(ssf*nl)



    F = np.zeros((Nx,Nx,nl, nt))


    x = np.arange(nx)
    X = np.linspace(0,nx-1,Nx)
    for t in range(nt):
        print(t)
        for l in range(nl):
            for j in range(ssf):
                fid = l + j*ssf + t*ssf*nl
                fname = p0_folder + 'p0_' + str(lambs[l]) + '/moby_' + str(fid) + '.mat'
                f = io.loadmat(fname)['p0'] 
                F[:,:,l,t] += np.mean(f,-1)/(ssf) 
    """try:
        print("Saving as 1 file")
        io.savemat('subsample_p0.mat',{'p0': F})"""
        
    """except:
        print("Saving as 3 files")"""
    io.savemat('subsample_p00.mat',{'p0': F[:,:,:,::3]}, do_compression = True)
    io.savemat('subsample_p01.mat',{'p0': F[:,:,:,1::3]}, do_compression = True)
    io.savemat('subsample_p02.mat',{'p0': F[:,:,:,2::3]}, do_compression = True)
        
        






    



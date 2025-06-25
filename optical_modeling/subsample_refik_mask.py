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
    
    anat_folder =  '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/dynamic_phantom_anatomy/'
    nx = 248
    zmin = 93
    zmax = 107


    lambs = [730,  750,  770,  780,  790,  800,  820,  840]
    nl = len(lambs)

    ssf = 1

    nT = 9600
    nt = nT//(ssf*nl)



    F1 = np.zeros((nx,nx, nt))
    F2 = np.zeros((nx,nx, nt))
    G = np.zeros((nx,nx,nt))


    for t in range(nt):
        print(t)
        for j in range(nl*ssf):
            fid = t*nl*ssf + j
            #fname = anat_folder + 'dynamic_phantom_anatomy_1.npy'
            fname = anat_folder + 'dynamic_phantom_anatomy_' + str(fid%199 + 1) + '.npy'
            """if fid%199 == 0:
                fname = anat_folder + 'dynamic_phantom_anatomy_' + str(fid%199 + 2) + '.npy'
            else:
                fname = anat_folder + 'dynamic_phantom_anatomy_' + str(fid%199 + 1) + '.npy'"""
            m = np.load(fname)[:,:,zmin:zmax]
            if t == 0:
                io.savemat('anat.mat',{'anat': m}, do_compression = True)
                assert False

            F1 [:,:,t] += np.sum(((m==97) + (m == 98)).astype(int),-1)
            F2 [:,:,t] += np.sum((m==98).astype(int),-1)
            G[:,:,t] += np.sum((m==96).astype(int),-1)
    F1 = (F1 > 0).astype(int)
    F2 = (F2 > 0).astype(int)
    G = (G > 0).astype(int)
    io.savemat('subsample_mask.mat',{'mask': F1, 'core_mask': F2}, do_compression = True)
    io.savemat('subsample_aif_mask.mat',{'aif_mask': G}, do_compression = True)
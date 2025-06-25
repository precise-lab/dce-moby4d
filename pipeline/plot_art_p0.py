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

    M = io.loadmat('../optical_modeling/subsample_aif_mask.mat')['aif_mask']
    M[:190,:,:] = 0
    M = np.transpose(M,(1,0,2))
    #plt.imshow(M[:,:,0])
    #plt.savefig("act_test.png")
    #assert False
    p0 = np.zeros((248,248, 8, 1200))
    p0[:,:,:,::3] = io.loadmat('../acoustic_modeling/subsample_p00.mat')['p0']
    p0[:,:,:,1::3] = io.loadmat('../acoustic_modeling/subsample_p01.mat')['p0']
    p0[:,:,:,2::3] = io.loadmat('../acoustic_modeling/subsample_p02.mat')['p0']

    aif_curves = np.sum( p0*M[:,:,None,:], (0,1))
    
    plt.clf()
    for j in range(8):
        plt.plot(aif_curves[j,:])
    plt.savefig('act_test.png')
    


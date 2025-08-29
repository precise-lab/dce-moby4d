import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine phantom elements', fromfile_prefix_chars='@')
    parser.add_argument('--zmin', default = 175, type = int)
    parser.add_argument('--zmax', default = 375, type = int)
    args = parser.parse_args()


    with open( '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/moby_anatomy/moby_act_1.bin', 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        phan = np.fromfile(file, dtype=np.float32)
    phan = phan.reshape((725, 248, 248))
    phan = np.transpose(phan, (2,1,0))

    with open('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyvein_r15.DAT', 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        vein = np.fromfile(file, dtype=np.float32)

    # Reshape the array to the desired dimensions (248, 248, 725)
    vein = vein.reshape((725, 248, 248))
    vein = np.transpose(vein, (2,1,0))

    with open('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyartery_r15.DAT', 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        artery = np.fromfile(file, dtype=np.float32)

    # Reshape the array to the desired dimensions (248, 248, 725)
    artery = artery.reshape((725, 248, 248))
    artery = np.transpose(artery, (2,1,0))

    tumor = loadmat('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_phantom.mat')
    tumor = tumor['tumor_phantom']

    phan[vein == 1] = 95
    phan[artery == 1] = 96
    phan[tumor == 1] = 97
    phan[tumor == 2] = 98

    zmin = args.zmin
    zmax = args.zmax

    np.save('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/phantom_anatomy_0.npy', phan[:,:, zmin:zmax])
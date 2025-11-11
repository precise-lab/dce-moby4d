import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

import h5py

import datetime

if __name__ == "__main__":

    input_shape = (725, 248, 248)

    moby_file = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/moby_anatomy/moby_act_1.bin'
    with open(moby_file , 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        phan = np.fromfile(file, dtype=np.float32)
    phan = phan.reshape(input_shape)
    phan = np.transpose(phan, (2,1,0))

    vein_file = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyvein_r15.DAT'
    with open(vein_file, 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        vein = np.fromfile(file, dtype=np.float32)

    # Reshape the array to the desired dimensions (248, 248, 725)
    vein = vein.reshape(input_shape)
    vein = np.transpose(vein, (2,1,0))

    artery_file = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyartery_r15.DAT'
    with open(artery_file, 'rb') as file:
        # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        artery = np.fromfile(file, dtype=np.float32)

    # Reshape the array to the desired dimensions (248, 248, 725)
    artery = artery.reshape(input_shape )
    artery = np.transpose(artery, (2,1,0))

    tumor_file = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_phantom.mat'
    tumor = loadmat(tumor_file)
    tumor = tumor['tumor_phantom']

    inside = (phan != 0)*(phan != 0) #remove background and skin
    phan[ (vein == 1) * inside ] = 95
    phan[ (artery == 1)* inside] = 96
    phan[ (tumor == 1)*inside  ] = 97
    phan[ (tumor == 2)*inside ] = 98

    phan = phan.astype(np.uint8)

    with h5py.File("phantom_anatomy_0.h5", "w") as fid:
        group = fid.create_group("Info")
        group.attrs["moby_file"]   = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/moby_anatomy/moby_act_1.bin'
        group.attrs["vein_file"]   = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyvein_r15.DAT'
        group.attrs["artery_file"] = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyartery_r15.DAT'
        group.attrs["tumor_file"]  = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_phantom.mat'
        group.attrs["created on"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        dset = fid.create_dataset("phantom", data=phan, compression="lzf")
        dset.attrs["spacing"] = np.array([0.15]*3, dtype=np.float32)
        dset.attrs["units"] = "mm"
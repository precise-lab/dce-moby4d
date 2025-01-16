import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

zmin = 175
zmax = 375


for i in range(1, 2):
    print(i)
    with open( '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/moby_anatomy/moby_act_{}.bin'.format(i), 'rb') as file:
    # Read the file as a float32 array (equivalent to 'float' in MATLAB)
        phan = np.fromfile(file, dtype=np.float32)

    # Reshape the array to the desired dimensions (248, 248, 725)
    phan = phan.reshape((725, 248, 248))
    phan = np.transpose(phan, (2,1,0))

    vein_moved = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/vein_moved/vein_moved_{}.npy'.format(i))
    artery_moved = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/artery_moved/artery_moved_{}.npy'.format(i))
    tumor_moved = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_moved/tumor_moved_{}.npy'.format(i))

    phan[vein_moved == 1] = 95
    phan[artery_moved == 1] = 96
    phan[tumor_moved == 1] = 97
    phan[tumor_moved == 2] = 98

    np.save('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/dynamic_phantom_anatomy/dynamic_phantom_anatomy_{}.npy'.format(i), phan[:,:, zmin:zmax])  
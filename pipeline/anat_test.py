import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io 

if __name__ == "__main__":
    vein = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/vein_moved/vein_moved_1.npy')
    io.savemat( 'test_vein.mat', {'vein':vein})

    artery = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/artery_moved/artery_moved_1.npy')
    io.savemat('test_artery.mat', {'artery':artery})

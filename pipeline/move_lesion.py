import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

tumor = loadmat('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_phantom.mat')
label_tumor = tumor['tumor_phantom']


for i in range(1, 200):
    print(i)
    tumor_moved = np.zeros((248, 248, 725))
    phan_motion = np.zeros((248, 248, 725))
    ones_array = np.ones((248, 248, 725))

    #frame1 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/motion_vectors_converted_to_numpy_arrays/frame1_{}.npy'.format(i))
    #frame2 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/motion_vectors_converted_to_numpy_arrays/frame{}.npy'.format(i))
    frame1 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame1_{}.npy'.format(i+1))
    frame2 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame{}.npy'.format(i+1))

    for j in range(frame1.shape[0]):
        if label_tumor[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] == 1.:
            tumor_moved[round(frame2[j, 0]), round(frame2[j, 1]), round(frame2[j, 2])] = 1
        if label_tumor[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] == 2.:
            tumor_moved[round(frame2[j, 0]), round(frame2[j, 1]), round(frame2[j, 2])] = 2
        phan_motion[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] = 1

    diff = ones_array - phan_motion
    tumor_moved = tumor_moved + diff*label_tumor

    #np.save('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/tumor_moved/tumor_moved_{}.npy'.format(i), tumor_moved)
    np.save('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_moved/tumor_moved_{}.npy'.format(i), tumor_moved)
    
import numpy as np
import matplotlib.pyplot as plt


# Open the .DAT file in binary mode
with open('onlyartery_r15.DAT', 'rb') as file:
    # Read the file as a float32 array (equivalent to 'float' in MATLAB)
    artery = np.fromfile(file, dtype=np.float32)

# Reshape the array to the desired dimensions (248, 248, 725)
artery = artery.reshape((725, 248, 248))
artery = np.transpose(artery, (2,1,0))
print(artery.shape)


for i in range(1, 200):
    print(i)
    artery_moved = np.zeros((248, 248, 725))

    # Some parts of the anatomical phantom do not have motion; therefore, they are not outputted in the motion vectors.
    # Nonetheless, there might be vessels in those parts of the phantom.
    # 'phan_motion' represents the regions where there is motion.
    phan_motion = np.zeros((248, 248, 725))

    #frame1 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/frame1_{}.npy'.format(i))
    #frame2 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/frame{}.npy'.format(i))
    frame1 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame1_{}.npy'.format(i+1))
    frame2 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame{}.npy'.format(i+1))

    for j in range(frame1.shape[0]):
        if artery[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] == 1:
            artery_moved[round(frame2[j, 0]), round(frame2[j, 1]), round(frame2[j, 2])] = 1
        phan_motion[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] = 1
    
    artery_masked = phan_motion*artery
    artery_diff = artery - artery_masked
    artery_moved = artery_moved + artery_diff

    #np.save('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/artery_moved_{}.npy'.format(i), artery_moved)
    np.save('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/artery_moved/artery_moved_{}.npy'.format(i), artery_moved)
    
    '''
    plt.figure()
    plt.imshow(np.squeeze(np.amax(artery_moved, axis=0)), vmin=0., vmax=1.)
    plt.savefig('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body/artery_moved_{}.png'.format(i))
    '''
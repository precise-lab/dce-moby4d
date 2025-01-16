import numpy as np
import matplotlib.pyplot as plt


# Open the .DAT file in binary mode
with open('onlyvein_r15.DAT', 'rb') as file:
    # Read the file as a float32 array (equivalent to 'float' in MATLAB)
    vein = np.fromfile(file, dtype=np.float32)

# Reshape the array to the desired dimensions (248, 248, 725)
vein = vein.reshape((725, 248, 248))
vein = np.transpose(vein, (2,1,0))
print(vein.shape)


for i in range(1, 200):
    print(i)
    vein_moved = np.zeros((248, 248, 725))

    # some part of the anatomical phantom doesn't have motion; therefore it is not outputed in the motion vectors
    # nontheless; there might be vessels in those parts of the phantom. Phan_motion is the region where there is motion
    phan_motion = np.zeros((248, 248, 725))

    #frame1 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/frame1_{}.npy'.format(i))
    #frame2 = np.load('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/frame{}.npy'.format(i))
    frame1 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame1_{}.npy'.format(i+1))
    frame2 = np.load('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame{}.npy'.format(i+1))

    for j in range(frame1.shape[0]):
        if vein[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] == 1:
            vein_moved[round(frame2[j, 0]), round(frame2[j, 1]), round(frame2[j, 2])] = 1
        phan_motion[round(frame1[j, 0]), round(frame1[j, 1]), round(frame1[j, 2])] = 1
    
    vein_masked = phan_motion*vein
    vein_diff = vein - vein_masked
    vein_moved = vein_moved + vein_diff

    #np.save('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body_anesthesia/vein_moved_{}.npy'.format(i), vein_moved)
    np.save('/workspace/shared_data/Moby_multi_wave/Refik_Mouse/vein_moved/vein_moved_{}.npy'.format(i), vein_moved)
    
    '''
    plt.figure()
    plt.imshow(np.squeeze(np.amax(vein_moved, axis=0)), vmin=0., vmax=1.)
    plt.savefig('/home/rcam2/Downloads/MOBY_phantom_generation/anatomical_structure_body/vein_moved_{}.png'.format(i))
    '''
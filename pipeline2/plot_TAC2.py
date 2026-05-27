import h5py 
import numpy as np
import matplotlib.pyplot as plt
import glob

from scipy.ndimage import center_of_mass



def averaged_rois(n_resp_phases):
    artery_volume = np.zeros(n_resp_phases)
    tumor_volume = np.zeros(n_resp_phases)
    tumor_core_volume = np.zeros(n_resp_phases)

    for resp_phase in np.arange(n_resp_phases):
        fname = f"indicator/indicator_{resp_phase:05d}.h5"
        with h5py.File(fname, "r") as fid:
            artery_tmp = fid["artery"][:,:,:]
            tumor_tmp = fid["tumor"][:,:,:]
            tumor_core_tmp = fid["tumor_core"][:,:,:]

        if resp_phase == 0:
            artery_roi = artery_tmp
            tumor_roi = tumor_tmp
            tumor_core_roi = tumor_core_tmp 
        else:
            artery_roi = artery_roi+artery_tmp
            tumor_roi = tumor_roi+tumor_tmp
            tumor_core_roi = tumor_core_roi+tumor_core_tmp 

        artery_volume[resp_phase] = np.sum(artery_tmp)
        tumor_volume[resp_phase] = np.sum(tumor_tmp)
        tumor_core_volume[resp_phase] = np.sum(tumor_core_tmp)

    plt.figure()
    plt.plot(artery_volume)
    plt.savefig("Arteryvolume.png")
    plt.figure()
    plt.plot(tumor_volume)
    plt.savefig("Tumorvolume.png")
    plt.figure()
    plt.plot(tumor_core_volume)
    plt.savefig("TumorCorevolume.png")

        

    artery_roi = artery_roi/float(n_resp_phases)
    tumor_roi = tumor_roi/float(n_resp_phases)
    tumor_core_roi = tumor_core_roi/float(n_resp_phases)

    delta = 40
    center = center_of_mass(tumor_core_roi)
    print(center)
    ranges_xyz = [ slice( [int(c-delta), int(c+delta)] ) for c in center[:] ]

    return ranges_xyz, artery_roi, tumor_roi, tumor_core_roi



if __name__=="__main__":
    nframes = 12*60*20+1
    tumor_tac = np.nan*np.ones(nframes)
    tumor_core_tac = np.nan*np.ones(nframes)
    aif       = np.nan*np.ones(nframes)
    time      = np.nan*np.ones(nframes)
    n_resp_phases = 200
    
    avg_ranges_xyz, avg_artery_roi, avg_tumor_roi, avg_tumor_core_roi = averaged_rois(n_resp_phases)

    for resp_phase in np.arange(n_resp_phases):
        fname = f"indicator/indicator_{resp_phase:05d}.h5"
        with h5py.File(fname, "r") as fid:
            artery_roi = fid["artery"][:,:,:]
            tumor_roi = fid["tumor"][:,:,:]
            tumor_core_roi = fid["tumor_core"][:,:,:]
        
        delta = 40
        center = center_of_mass(tumor_core_roi)
        ranges_xyz = [ slice( [int(c-delta), int(c+delta)] ) for c in center[:] ]
        artery_roi = artery_roi[*ranges_xyz]
        tumor_roi = tumor_roi[*ranges_xyz]
        tumor_core_roi = tumor_core_roi[*ranges_xyz]
                                               

        fname = f"p0/p0_{resp_phase:05d}.h5"
        with h5py.File(fname, "r") as fid:
             p0 = fid["p0"][ranges_xyz[0],
                            ranges_xyz[1],
                            ranges_xyz[2],
                            :]
             time_ll = fid["time"][:]
             for local_index in range(time_ll.shape[0]):
                time_index = local_index*n_resp_phases + resp_phase
                time[time_index] = time_ll[local_index]
                aif[time_index] = np.sum(p0*artery_roi) / np.sum(artery_roi) #np.percentile(p0[artery_roi>0.5,local_index], 99) #np.sum(p0*avg_artery_roi) / np.sum(avg_artery_roi)
                tumor_tac[time_index] = np.sum(p0*tumor_roi) / np.sum(tumor_roi) #np.percentile(p0[tumor_roi>0.5, local_index], 99)#np.sum(p0*avg_tumor_roi) / np.sum(avg_tumor_roi)
                tumor_core_tac[time_index] = np.sum(p0*tumor_core_roi) / np.sum(tumor_core_roi)#np.percentile(p0[tumor_core_roi>0.5, local_index], 99)#np.sum(p0*avg_tumor_core_roi) / np.sum(avg_tumor_core_roi)

                print(f"{time[time_index]:0.2f}, {aif[time_index]:0.5f}, {tumor_tac[time_index]:0.5f}, {tumor_core_tac[time_index]:0.5f}")

    is_nan = np.isnan(time)
    time = time[~is_nan]
    aif = aif[~is_nan]
    tumor_tac = tumor_tac[~is_nan]
    tumor_core_tac = tumor_core_tac[~is_nan]

    plt.figure()
    plt.plot(time, aif)
    plt.savefig("Aif.png")

    plt.figure()
    plt.plot(time, tumor_tac)
    plt.savefig("Tumor_TAC.png")

    plt.figure()
    plt.plot(time, tumor_tac)
    plt.savefig("Tumor_CORE_TAC.png")



import h5py 
import numpy as np
import matplotlib.pyplot as plt
import glob

from scipy.ndimage import center_of_mass



def averaged_rois(n_resp_phases, phases_computed):
    artery_volume = np.zeros(phases_computed)
    tumor_volume = np.zeros(phases_computed)
    tumor_core_volume = np.zeros(phases_computed)

    for resp_phase in np.arange(phases_computed):
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

        

    artery_roi = artery_roi/float(phases_computed)
    tumor_roi = tumor_roi/float(phases_computed)
    tumor_core_roi = tumor_core_roi/float(phases_computed)

    center = center_of_mass(tumor_core_roi)
    print(center)

    radius = 50
    x, y, z = np.indices(artery_roi.shape)
    # Euclidean distance from center
    dist = np.sqrt(
            (x - center[0])**2 +
            (y - center[1])**2 +
            (z - center[2])**2
        )

    # set points at distance == radius
    arr = np.zeros_like(artery_roi)
    arr[np.isclose(dist, radius, atol=0.5)] = 1

    return arr, artery_roi, tumor_roi, tumor_core_roi



if __name__=="__main__":
    nframes = 12*60*20
    tumor_tac = np.nan*np.ones(nframes)
    tumor_core_tac = np.nan*np.ones(nframes)
    aif       = np.nan*np.ones(nframes)
    time      = np.nan*np.ones(nframes)
    n_resp_phases = 200
    phases_computed = 71
    
    avg_roi, avg_artery_roi, avg_tumor_roi, avg_tumor_core_roi = averaged_rois(n_resp_phases, phases_computed)
    avg_artery_roi = avg_roi*avg_artery_roi

    for resp_phase in np.arange(phases_computed):
        fname = f"indicator/indicator_{resp_phase:05d}.h5"
        with h5py.File(fname, "r") as fid:
            artery_roi = fid["artery"][:,:,:]
            tumor_roi = fid["tumor"][:,:,:]
            tumor_core_roi = fid["tumor_core"][:,:,:]

        artery_roi = artery_roi * avg_roi

        for time_index in np.arange(0+resp_phase, nframes, n_resp_phases):
            fname = f"p0/p0_{time_index:05d}.h5"
            with h5py.File(fname, "r") as fid:
                p0_set = fid["p0"]
                p0 = p0_set[:,:,:]

                time[time_index] = p0_set.attrs["time"]
                aif[time_index] = np.percentile(p0[artery_roi>0.5], 99) #np.sum(p0*avg_artery_roi) / np.sum(avg_artery_roi)
                tumor_tac[time_index] =  np.percentile(p0[tumor_roi>0.5], 99)#np.sum(p0*avg_tumor_roi) / np.sum(avg_tumor_roi)
                tumor_core_tac[time_index] = np.percentile(p0[tumor_core_roi>0.5], 99)#np.sum(p0*avg_tumor_core_roi) / np.sum(avg_tumor_core_roi)

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



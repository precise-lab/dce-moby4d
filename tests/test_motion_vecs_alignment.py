import h5py
import argparse
import numpy as np
from scipy.ndimage import shift, binary_dilation

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='test_alignement', fromfile_prefix_chars='@')
    parser.add_argument('-p', '--phantom', default = "phantom_anatomy_0.h5")
    parser.add_argument('-d', '--disp', default = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors.h5')

    
    args = parser.parse_args()

    with h5py.File(args.phantom, "r") as fid:
        data_set = fid["/phantom"]
        h = data_set.attrs["spacing"]
        data = data_set[:,:,:].astype(np.int32)

    with h5py.File(args.disp) as fid:
            iset = fid["interpolated"]
            imposed_array = (iset[:,:,:,0] == 0)
            for iframe in range(1, iset.shape[-1]):
                imposed_array = np.logical_or( iset[:,:,:,iframe] == 0, imposed_array)

    imposed_array_dilated = binary_dilation(imposed_array, iterations=2)
    imposed_array_shifted = shift(imposed_array, shift=(1,1,1), cval = False)
    imposed_organs = np.unique( data[imposed_array_dilated].flatten() )

    for organ in imposed_organs[:]:
         volume = np.sum(data == organ )
         moved_vol_dilated = np.sum(data[imposed_array_dilated] == organ )
         moved_vol = np.sum(data[imposed_array] == organ )
         moved_vol_s = np.sum(data[imposed_array_shifted] == organ )

         print(f"{organ}  {volume} {moved_vol_dilated} {moved_vol} {moved_vol_s}")

    
import numpy as np


def run(nframes, N, h, fname_frame1, fname_cframe, fname_out):
	for iframe in range(1, nframes):
		print(iframe)
		# Load the frames
		frame_0 = np.load(fname_frame1.format(iframe+1))
		frame_current = np.load( fname_cframe.format(iframe+1))

		# Calculate the motion vector differences using NumPy's vectorized operations
		motion_vectors = frame_current - frame_0
		motion_vectors[:,0] *= h[0]
		motion_vectors[:,1] *= h[1]
		motion_vectors[:,2] *= h[2]

		motion = np.zeros(N+[3]) #Create 4D tensor
		for c in np.arange(frame_0.shape[0]):
			i = int(np.floor( frame_0[c,0] - 0.1 ))
			j = int(np.floor( frame_0[c,1] - 0.1 ))
			k = int(np.floor( frame_0[c,2] - 0.1 ))
			motion[i,j,k,:] = motion_vectors[c, :] 

		# Save the result as a NumPy array
		np.savez_compressed(fname_out.format(iframe), motion=motion)

if __name__ == "__main__":
	nframes = 200
	N = [248, 248, 725]
	h = [0.15, 0.15, 0.15]
	fname_frame1 = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame1_{0}.npy'
	fname_cframe = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/frame{0}.npy'
	fname_out    = '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_array/motion{0}.npz'

	run(nframes, N, h,fname_frame1, fname_cframe, fname_out)


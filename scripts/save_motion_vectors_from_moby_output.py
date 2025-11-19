import csv
import numpy as np
import h5py
import os
import datetime


def process_frame(fname, iframe, N, h):
	print("Processing frame: ", iframe)
	with open(fname, 'r') as file:
		data = file.read()

	# Splitting the data into lines
	lines = data.strip().split('\n')
	assert lines[0].strip() == "Surface Vectors From Frame 1 to Frame {}".format(iframe+1)

	frame_ref_pos = []
	frame_cur_pos = []
	interpolated = []
	for line in lines[1:]:
		row = line.split()
		assert row[1].strip() == "frame1"
		assert row[5].strip() == "frame{}".format(iframe+1)
		interpolated.append(row[0] == "interpolated_vector")
		frame_ref_pos.append([float(d) for d in row[2:5]])
		frame_cur_pos.append([float(d) for d in row[6:9]])

	frame_ref_pos = np.array(frame_ref_pos, dtype=np.float64)
	frame_cur_pos = np.array(frame_cur_pos, dtype=np.float64)
	interpolated_v  = np.array(interpolated, dtype=np.int8)

	print("Read frame positions of size:", frame_ref_pos.shape)
	print("Smallest index: ", np.min(frame_ref_pos, axis=0))
	print("Largest index:  ", np.max(frame_ref_pos, axis=0))

	motion_vectors = frame_cur_pos - frame_ref_pos
	motion       = np.zeros(N+[3], dtype=np.float32) #Create 4D tensor
	interpolated = -1 * np.ones(N, dtype=np.int8)
	for c in np.arange(frame_ref_pos.shape[0]):
		i = int(np.floor( frame_ref_pos[c,0] - 0.1 ))
		j = int(np.floor( frame_ref_pos[c,1] - 0.1 ))
		k = int(np.floor( frame_ref_pos[c,2] - 0.1 ))
		motion[i,j,k,:] = motion_vectors[c, :]*h
		interpolated[i,j,k] = interpolated_v[c]

	return motion, interpolated

if __name__=="__main__":
	fname = "motion_vectors_text_moby_output/moby_vec_frame1_to_frame{}.txt"
	output = "motion_vectors.h5"
	N = [248, 248, 725]
	h = np.array([0.15, 0.15, 0.15])
	dt = 0.05
	nframes = 200 

	m_size = N + [3, nframes]
	m_chunks = N + [3, 1]
	i_size = N + [nframes]
	i_chunks = N + [1]
	with h5py.File(output, "w") as fid:
		group = fid.create_group("Info")
		group.attrs["created from"] = os.path.abspath("motion_vectors_text_moby_output/moby_vec_frame1_to_frame{}.txt")
		group.attrs["created on"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		mset = fid.create_dataset("motion", m_size, chunks=tuple(m_chunks), compression="lzf", dtype=np.float32)
		iset = fid.create_dataset("interpolated", i_size, chunks=tuple(i_chunks), compression="lzf", dtype=np.int8)
		mset.attrs["dt"] = 0.05
		mset.attrs["units"] = "mm"
		mset.attrs["created from"] = os.path.abspath("motion_vectors_text_moby_output/moby_vec_frame1_to_frame{}.txt")
		mset.attrs["spacing"] = h
		mset[:,:,:,:,0] = np.zeros(m_chunks[0:-1], dtype=np.float32)
		iset[:,:,:,0] = -np.ones(i_chunks[0:-1], dtype=np.int8)
		for iframe in range(1,nframes):
			motion, interpolated = process_frame(fname.format(iframe+1), iframe, N, h)
			mset[:,:,:,:,iframe] = motion
			iset[:,:,:,iframe] = interpolated
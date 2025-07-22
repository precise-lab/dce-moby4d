export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

mpiexec -n 4 python -u combined_optical_forward.py -w 840 --wb 7 -p /workspace/shared_data/Moby_multi_wave/p0_840/ --pzmin 93 --pzmax 107
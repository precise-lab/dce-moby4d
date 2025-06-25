export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

mpiexec -n 8 python -u combined_optical_forward.py -w 820 --wb 6 -p /workspace/shared_data/Moby_multi_wave/p0_820/
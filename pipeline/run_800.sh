export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

mpiexec -n 4 python -u combined_optical_forward.py -w 800 --wb 5 -p /workspace/shared_data/Moby_multi_wave/p0_800/ -s /workspace/shared_data/Moby_multi_wave/so2/ --pzmin 93 --pzmax 107
#mpiexec -n 8 python -u form_mu_a.py -w 800 --wb 5 -p /workspace/shared_data/Moby_multi_wave/mu_a_800/
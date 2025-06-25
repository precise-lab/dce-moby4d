export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

mpiexec -n 8 python -u combined_optical_forward.py -w 790 --wb 4 -p /workspace/shared_data/Moby_multi_wave/p0_790/ -s /workspace/shared_data/Moby_multi_wave/so2/
#mpiexec -n 8 python -u form_mu_a.py -w 790 --wb 4 -p /workspace/shared_data/Moby_multi_wave/mu_a_790/
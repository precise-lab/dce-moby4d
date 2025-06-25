export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

mpiexec -n 8 python -u combined_optical_forward.py -w 770 --wb 2 -p /workspace/shared_data/Moby_multi_wave/p0_770/ 
#mpiexec -n 8 python -u form_mu_a.py -w 750 --wb 2 -p /workspace/shared_data/Moby_multi_wave/mu_a_770/ 
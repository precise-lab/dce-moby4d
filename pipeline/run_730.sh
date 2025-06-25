export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1


mpiexec -n 4 python -u combined_optical_forward.py -w 730 --wb 0 -p /workspace/shared_data/Moby_multi_wave/p0_730/
#mpiexec -n 8 python -u form_mu_a.py -w 730 --wb 0 -p /workspace/shared_data/Moby_multi_wave/mu_a_730/
export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

python compute_displacement_field.py
python combine_phantom_elements.py --zmin 175 --zmax 375
python generate_mesh_structured.py 
python move_mesh.py
mpiexec -n 4 python -u compute_so2.py 
mpiexec -n 4 python -u compute_mask.py --pzmin 93 --pzmax 107

./run_730.sh
./run_750.sh
./run_770.sh
./run_790.sh
./run_800.sh
./run_820.sh
./run_840.sh
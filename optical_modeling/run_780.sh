export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1
export OMP_NUM_THREADS=1

python -u solve_pa.py -l 780 -s 0 -e 600
python -u solve_pa.py -l 780 -s 600 -e 1200
python -u solve_pa.py -l 780 -s 1200 -e 1800
python -u solve_pa.py -l 780 -s 1800 -e 2400
python -u solve_pa.py -l 780 -s 2400 -e 3000
python -u solve_pa.py -l 780 -s 3000 -e 3600
python -u solve_pa.py -l 780 -s 3600 -e 4200
python -u solve_pa.py -l 780 -s 4200 -e 4800
python -u solve_pa.py -l 780 -s 4800 -e 5400
python -u solve_pa.py -l 780 -s 5400 -e 6000
python -u solve_pa.py -l 780 -s 6000 -e 6600
python -u solve_pa.py -l 780 -s 6600 -e 7200
python -u solve_pa.py -l 780 -s 7200 -e 7800
python -u solve_pa.py -l 780 -s 7800 -e 8400
python -u solve_pa.py -l 780 -s 8400 -e 9000
python -u solve_pa.py -l 780 -s 9000 -e 9600
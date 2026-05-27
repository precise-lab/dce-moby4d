for start in $(seq 0 10 199); do
    CMD="sbatch --export=ALL,start=$start job.slurm"
    echo $CMD
done

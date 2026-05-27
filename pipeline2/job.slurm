#!/bin/bash
#SBATCH -J moby           # Job name
#SBATCH -o moby.%j.out        # Output file (%j = job ID)
#SBATCH -t 05:00:00              # Wall time (1 hour)
#SBATCH -N 2                     # Number of nodes
#SBATCH -n 128                     # Total tasks (4 per node)
#SBATCH -p gg                    # GPU partition (modify as needed)
#SBATCH -A ASC25089          # Project allocation

source $WORK/fenics-env/bin/activate

which python

# Create directory if it does not exist
if [ ! -d p0 ]; then
    mkdir -p p0
fi

# Create directory if it does not exist
if [ ! -d indicator ]; then
    mkdir -p indicator
fi

# Create directory if it does not exist
if [ ! -d logs ]; then
    mkdir -p logs
fi

printf -v padded_start "%04d" "$start"
rm -rf rerun_$padded_start

export UCX_TLS=rc,sm,self
export UCX_NET_DEVICES=mlx5_0:1
export UCX_RC_TIMEOUT=10s
export UCX_RC_RETRY_COUNT=20

end=$((start + 9))
date
for P in $(seq $start $end); do
    printf -v padded_P "%04d" "$P"
    CMD="ibrun python -u run_optical_model2.py -p $P >& logs/run_optical_model2_$padded_P.log"

    # Print the command before executing
    echo "$CMD"

    # Execute the command
    ibrun python -u run_optical_model2.py -p $P >& logs/run_optical_model2_$padded_P.log

    if [ $? -ne 0 ]; then
	   echo  "Failed $P"
	   echo "$CMD" >> rerun_$padded_start
    fi
done
date

#!/bin/bash
#SBATCH --job-name=ferroSim
#SBATCH --output=sp-db1.out
#SBATCH --error=sp-db1.err
#SBATCH --partition=share
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G

totalExecCores=$(($SLURM_NTASKS_PER_NODE*$SLURM_CPUS_PER_TASK*$SLURM_JOB_NUM_NODES))
execMemory=$(($SLURM_MEM_PER_NODE/$SLURM_NTASKS_PER_NODE))M
echo Total Cores: $totalExecCores
echo Total Memory: $execMemory
echo CPUS per Task: $SLURM_CPUS_PER_TASK
echo Number of Nodes: $SLURM_JOB_NUM_NODES

export TMPDIR="$HOME/tmp"
export SPARK_LOCAL_DIRS=$(mktemp -d)
export SPARK_WORKER_DIR=$(mktemp -d)
export SPARK_LOG_DIR="$TMPDIR/spark-logs"
mkdir -p $SPARK_LOG_DIR

MASTER=$(hostname)
echo "MASTER ADDRESS: $MASTER:7077"
srun start-cluster.sh $MASTER &

sleep 30

#zip up folder for dependencies
rm dependencies.zip
zip -r dependencies.zip . &>/dev/null
echo "Zipped dependencies"

echo "About to submit job"
spark-submit \
	--class spark-python-ferro \
	--master spark://$MASTER:7077 \
        --num-executors $SLURM_JOB_NUM_NODES \
        --executor-cores $SLURM_CPUS_PER_TASK \
        --executor-memory $execMemory \
        --py-files dependencies.zip \
	ferroScript.py

tail -f /dev/null
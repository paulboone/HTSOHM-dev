#!/bin/bash
#
# Example shell script for running job that runs off the Wilmerlab subjobserver.
# $Revision: 1.0 $
# $Date:  2016-11-10 $
# $Author: akaija $

#PBS -j oe
#PBS -N htsohm
#PBS -q shared
#PBS -l nodes=1:ppn=20
#PBS -l walltime=10:00:00
#PBS -l mem=2GB
#PBS -S /bin/bash

# accepts a parameter stay_alive if you don't want the worker to exit immediately after all jobs
# are complete
# use like `qsub -v stay_alive=1`

run_id = ${input}

echo JOB_ID: $PBS_JOBID JOB_NAME: $PBS_JOBNAME HOSTNAME: $PBS_O_HOST
echo start_time: `date`

# dependencies
module purge
module load python/3.5.1
module load postgresql/9.5.2
source ~/venv/htsohm/bin/activate

cd $PBS_O_WORKDIR
for ((i = 0; i < $PBS_NUM_PPN; i++))
do
    ./hts.py launch_worker ${run_id} >> ${run_id}/output_${PBS_O_HOST}_$$_$i.log 2>&1 &
done

wait

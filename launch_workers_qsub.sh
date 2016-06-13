#!/bin/bash
#
# Example shell script for running job that runs off the Wilmerlab subjobserver.
# $Revision: 1.0 $
# $Date:  2016-03-21 $
# $Author: paulboone $

#PBS -j oe
#PBS -N htsohm
#PBS -q test
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:15:00
#PBS -l mem=1GB
#PBS -S /bin/bash

# accepts a parameter stay_alive if you don't want the worker to exit immediately after all jobs
# are complete
# use like `qsub -v stay_alive=1`

echo JOB_ID: $PBS_JOBID JOB_NAME: $PBS_JOBNAME HOSTNAME: $PBS_O_HOST
echo start_time: `date`

# dependencies
module purge
module load python/3.5.1
module load postgresql/9.5.2
source ~/venv/htsohm/bin/activate

cd $PBS_O_WORKDIR
sjs_launch_workers.sh $PBS_NUM_PPN $stay_alive

# workaround for .out / .err files not always being copied back to $PBS_O_WORKDIR
cp /var/spool/torque/spool/$PBS_JOBID.OU $PBS_O_WORKDIR/$PBS_JOBID$(hostname)_$$.out
cp /var/spool/torque/spool/$PBS_JOBID.ER $PBS_O_WORKDIR/$PBS_JOBID$(hostname)_$$.err

exit

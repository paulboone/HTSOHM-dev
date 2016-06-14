#!/usr/bin/env python3
"""
This script supports reentry, i.e. if it gets killed for any reason, you should be able to rerun it
to resume where you left off.

"""

from time import sleep
import os
import sys
from datetime import datetime

import sjs
from rq.queue import get_failed_queue
from rq.job import Job
from rq.registry import StartedJobRegistry
import yaml

from htsohm.htsohm_hpc import manage_run, start_run
from htsohm.runDB_declarative import Material, session

MANAGE_RUN_TTL = 60*60*24*30 #30 days
RESTART_FILENAME = 'restart.yaml'

sjs.load(os.path.join("settings","sjs.yaml"))
redis_conn = sjs.get_redis_conn()
sjs_config = sjs.get_sjs_config()

def requeue_failed_jobs():
    failed_queue = get_failed_queue(connection=redis_conn)

    jobs_to_requeue = failed_queue.get_job_ids()
    for job_id in jobs_to_requeue:
        failed_queue.requeue(job_id)

    return jobs_to_requeue

def jobs_running():
    registry = StartedJobRegistry(name=sjs_config['queue'], connection=redis_conn)
    return registry.get_job_ids()

def jobs_queued():
    jobs_queue = sjs.get_job_queue()
    return jobs_queue.get_job_ids()

def timestamp():
    return datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

# def last_generation():
#     return session.query(func.max(Material.generation)).filter(Material.run_id == run_id).one()[0]

def unfinished_materials(generation):
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.generation == generation,
        Material.data_complete == False
    ).all()

def manage_run_results(job_id):
    return Job.fetch(job_id, connection=redis_conn).result

if len(sys.argv) > 1:
    # we were passed new parameters, so we're starting a new run!
    prior_manage_job_id = None
    generation = 0
    run_id = start_run(
        int(sys.argv[1]),
        int(sys.argv[2]),
        float(sys.argv[3]),
        int(sys.argv[4])
    )
    print("[%s] starting new run with passed parameters" % timestamp())
else:
    # no parameters were passed, so we'll attempt to restart a previous run
    with open(RESTART_FILENAME, 'r') as yaml_file:
        restart_settings = yaml.load(yaml_file)

    prior_manage_job_id = restart_settings['prior_manage_job_id']
    generation          = restart_settings['generation']
    run_id              = restart_settings['run_id']
    print("[%s] resuming run with parameters in restart_settings file: %s"
        % (timestamp(), restart_settings))

if __name__ == "__main__":
    first_time_through = True

    while (True):

        if not first_time_through:
            sleep(60)
        else:
            first_time_through = False

        if prior_manage_job_id:
            results = manage_run_results(prior_manage_job_id)
            if results is None:
                # either hasn't finished or failed
                # in either case, if everything else clears, we'll try again at the end of the loop
                pass
            elif results == -1:
                print("The last manage_run job with id %s indicated the process is complete!")
                os.remove('restart.yaml')
                sys.exit(0)
            elif results > 0:
                generation = results

        if len(jobs_queued()) > 0:
            print("jobs still queued...")
            continue

        if len(jobs_running()) > 0:
            print("jobs still running...")
            continue

        if len(requeue_failed_jobs()) > 0:
            print("jobs just requeued...")
            continue

        # final check to make sure the results look like we expect in the DB:
        if len(unfinished_materials(generation)) > 0:
            raise SystemException("System looks like it is running normally, but we are missing results in the database!")

        # there are no failed jobs to requeue, no jobs in the queue and no jobs running, so we can
        # move on to the next generation
        job_queue = sjs.get_job_queue()
        prior_manage_job_id = job_queue.enqueue(
            manage_run,
            run_id, generation,
            result_ttl=MANAGE_RUN_TTL
        ).id

        print("[%s] queued manage_job with id = %s" % (timestamp(), prior_manage_job_id))

        # save restart file locally, for use by the user in case this process is killed
        with open(RESTART_FILENAME, 'w') as f:
            f.write(yaml.dump({
                'prior_manage_job_id': prior_manage_job_id,
                'generation': generation,
                'run_id': run_id
            }))

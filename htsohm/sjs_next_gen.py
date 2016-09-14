#!/usr/bin/env python3

from time import sleep
import os
import sys
from datetime import datetime

import sjs
from rq.queue import get_failed_queue
from rq.job import Job
from rq.registry import StartedJobRegistry
import yaml

from htsohm.htsohm_hpc import manage_run
from htsohm.runDB_declarative import Material, session

MANAGE_RUN_TTL = 60*60*24*30 #30 days
RESTART_FILENAME = 'restart.yaml'

sjs.load(os.path.join("settings","sjs.yaml"))
redis_conn = sjs.get_redis_conn()
sjs_config = sjs.get_sjs_config()

def manage_run_results(job_id):
    return Job.fetch(job_id, connection=redis_conn).result


with open(RESTART_FILENAME, 'r') as yaml_file:
    restart_settings = yaml.load(yaml_file)

prior_manage_job_id = restart_settings['prior_manage_job_id']
generation          = restart_settings['generation']
run_id              = restart_settings['run_id']

## look at results of last manage run job and adjust generation if necessary

if prior_manage_job_id:
    results = manage_run_results(prior_manage_job_id)
    if results is None:
        # either hasn't finished or failed
        # in either case, if everything else clears, we'll try again at the end of the loop
        pass
    elif results == -1:
        print("The last manage_run job with id %s indicated the process is complete!")
        os.remove('restart.yaml')
        sys.exit(64)
    elif results > 0:
        generation = results

## queue next generation
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

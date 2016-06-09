from time import sleep
import os

import yaml
from rq.registry import StartedJobRegistry
from sqlalchemy import func
import sjs

from htsohm.utilities import write_config_file, read_config_file
from htsohm.htsohm import seed_generation, next_generation, hpc_job_run_all_simulations
from htsohm.runDB_declarative import Material, session

sjs.load(os.path.join("settings","sjs.yaml"))
htsohm_queue = sjs.get_job_queue()
redis_conn = sjs.get_redis_conn()
sjs_config = sjs.get_sjs_config()


def start_run(
        children_per_generation,    # number of materials per generation
        number_of_atomtypes,        # number of atom-types per material
        strength_0,                 # intial strength parameter
        number_of_bins,             # number of bins for analysis
        max_generations=20):        # maximum number of generations

    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations)['run-id']

    print("starting run with run_id: %s" % run_id)
    job = htsohm_queue.enqueue(manage_run,run_id, timeout=sjs_config['max_seconds_per_job'])

def continue_run(run_id):
    htsohm_queue.enqueue(manage_run, run_id)

def manage_run(run_id):
    print("======================================================================================")
    print("== manage_run")
    print("======================================================================================")

    config = read_config_file(run_id)
    last_generation = session.query(func.max(Material.generation)).filter(Material.run_id == run_id).one()[0]

    if last_generation is None:
        ## this is the first time we've run!
        generation = 0
    else:
        # block until all jobs are complete
        # the only started job in the queue should be THIS job!
        registry = StartedJobRegistry(name='htsohm_queue', connection=redis_conn)
        running_job_ids = registry.get_job_ids()
        while len(running_job_ids) > 1:
            print("still %s unfinished jobs in the htsohm queue; sleeping 60s..." % (len(running_job_ids) - 1))
            sleep(60)
            running_job_ids = registry.get_job_ids()

        #
        # if we have any materials for this generation that are still not populated with data
        # we need to attempt to rerun them.
        unfinished_materials = session.query(Material).filter(Material.run_id == run_id,
                                                              Material.generation == last_generation,
                                                              Material.data_complete == False).all()

        if len(unfinished_materials) > 0:
            print("Found %s unfinished materials from the last generation!" % len(unfinished_materials))
            print("Requeueing them all: %s" % unfinished_materials)

            for material in unfinished_materials:
                htsohm_queue.enqueue(hpc_job_run_all_simulations, material.id,
                                     timeout=sjs_config['max_seconds_per_job'])

            htsohm_queue.enqueue(manage_run, run_id)
            return

        print ("No unfinished materials from the last_generation! Moving on to the next generation.")
        generation = last_generation + 1


    print("======================================================================================")
    print("== manage_run, generation = %s" % generation)
    print("======================================================================================")


    #
    # prepare, mutate, and queue up the next generation
    if generation > config['max-number-of-generations']:
        print("max generations exceeded; we're done here!")
        return
    elif generation == 0:
        # SEED GENERATION
        seed_generation(run_id,
                        config['children-per-generation'], config['number-of-atom-types'],
                        queue=htsohm_queue)
    elif generation >= 1:
        # FIRST GENERATION, AND ON...
        next_generation(run_id,
                        config['children-per-generation'], generation,
                        queue=htsohm_queue)

    # COMMENTING THIS OUT TO MAKE THIS MANUAL FOR NOW, PRIOR TO REFACTOR
    #
    # seed_generation or next_generation will have enqueued all the materials jobs, but now we
    # enqueue another manage_run job to run after they are all completed.
    # htsohm_queue.enqueue(manage_run, run_id)

    print("======================================================================================")
    print("== END manage_run, generation = %s" % generation)
    print("======================================================================================")

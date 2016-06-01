#!/usr/bin/env python3

import yaml
from redis import Redis
from rq import Worker, Queue

with open('redis.yaml', 'r') as yaml_file:
    hpc_config = yaml.load(yaml_file)

redis_conn = Redis(**hpc_config['redis'])

def job_string(j):
    if j is None:
        return "-"
    return "  %s (%s): %s" % (j.get_id(), j.get_status(), j.get_call_string())


print("## WORKERS ##")
ws = Worker.all(connection=redis_conn)
for w in sorted(ws, key=lambda x: x.name):
    print("worker %s: %s" % (w.name, job_string(w.get_current_job())))

print("")
print("## QUEUES ##")
qs = Queue.all(connection=redis_conn)
for q in sorted(qs, key=lambda x: x.name):
    print("queue %s:" % q.name)
    for j in q.get_jobs():
        print(job_string(j))

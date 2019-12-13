#!/usr/bin/env python3

from glob import glob
import click

from htsohm import htsohm_run

@click.command()
@click.argument('config_path', type=click.Path())
@click.option('-r', '--restart',  type=int, default=-1, help="generation to restart at")
@click.option('-f', '--override-database-errors',  is_flag=True, default=False)
@click.option('-n', '--num-processes', type=int, default=1)
def cmdline(config_path, restart, override_database_errors, num_processes):
    if config_path == "auto":
        yaml_files = glob("*.yaml")
        if len(yaml_files) == 0:
            print("ERROR: no YAML file in the current directory; `dps auto` cannot be used.")
            return
        config_path = yaml_files[0]
        if len(yaml_files) > 1:
            print("WARNING: more than one YAML file found in this directory. Using first one: %s" % config_path)

    htsohm_run(config_path, restart, override_database_errors, num_processes)

if __name__ == '__main__':
    cmdline()

#!/usr/bin/env python3

from glob import glob
import click

from htsohm.htsohm_serial import serial_runloop

@click.command()
@click.argument('config_path', type=click.Path())
def cmdline(config_path):
    if config_path == "auto":
        yaml_files = glob("*.yaml")
        if len(yaml_files) == 0:
            print("ERROR: no YAML file in the current directory; `dps auto` cannot be used.")
            return
        config_path = yaml_files[0]
        if len(yaml_files) > 1:
            print("WARNING: more than one YAML file found in this directory. Using first one: %s" % config_path)

    serial_runloop(config_path)

if __name__ == '__main__':
    cmdline()

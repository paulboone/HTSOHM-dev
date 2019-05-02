#!/usr/bin/env python3

import click

from htsohm.htsohm_serial import serial_runloop

@click.command()
@click.argument('config_path',type=click.Path())
def cmdline(config_path):
    serial_runloop(config_path)

if __name__ == '__main__':
    cmdline()

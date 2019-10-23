"""
simple io.StringIO based log intended to be used with different workers so that work is output
in full blocks, where a block is a complete unit of work for the worker.
"""

import io

__buffered_log__ = ""

def init_slog():
    global __buffered_log__
    __buffered_log__ = io.StringIO("")
    return

def slog(*args):
    print(*args, file=__buffered_log__)

def get_slog_file():
    return __buffered_log__

def get_slog():
    return __buffered_log__.getvalue()

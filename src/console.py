"""
This file implements a simple console log
"""

from .container import Container

LOG = ""
SUBLOG = ""

HAS_SUBLOG = False


def log_print(msg):
    global LOG, SUBLOG, HAS_SUBLOG
    print(msg)
    LOG += f"\n{msg}"
    if HAS_SUBLOG:
        SUBLOG += f"\n{msg}"


def start_sublog():
    global SUBLOG, HAS_SUBLOG
    SUBLOG = ""
    HAS_SUBLOG = True


def export_sublog(container: Container):
    global SUBLOG, HAS_SUBLOG
    log_path = container.log_dir.log_path
    log_file = open(log_path, 'w+')
    log_file.write(SUBLOG)
    log_file.close()
    SUBLOG = None
    HAS_SUBLOG = False

def get_sublog():
    global SUBLOG
    return SUBLOG

def show_log():
    global LOG
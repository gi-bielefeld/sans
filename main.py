"""
This file implements the driver code for the testing software
and defines the software entry point
@Author Fabian Kolesch
"""

import argparse

from src import console
from src.config import Config
from src import script_parser
from src import time

from src.dir import *
from src import pipeline_factory

from src import seq_make


if __name__ == '__main__':
    log_print("--- RUNNING SANS_CI ---")
    log_print("\n[MAIN] Bootstrapping")
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", help="The target run script")
    parser.add_argument("--dry", help="Parse the run script without execution", default=False, action="store_true")
    args = parser.parse_args()

    # Load Config file
    Config.init()

    # Get the run script
    run_script = args.run
    # Create run containers
    containers = script_parser.parse_run_script(run_script)
    # If this was a test run, finish
    if args.dry:
        exit(0)

    log_print("\n[MAIN] Launching Tests")
    stamp = time.get_stamp()
    for container in containers:
        if container.type == "truth":
            log_print(f"Truth testing is currently disabled, due to reult inconsistencies ... skipping container {container.name}")
            continue
        console.start_sublog()
        log_print(f"\n\n[PROCESSING CONTAINER] {container.name} at STAMP: {stamp}")
        log_print(f"{container.get_info()}")
        seq_dir = Seq(container)
        run_dir = RunData(container)
        log_dir = TestLog(container)

        # Data creation
        if container.data.source == "synthetic":
            seq_make.make_synthetic_truth(container)
        
        elif container.data.source == "comp":
            seq_make.make_comp_data(container)

        pipeline_fac = None
        if container.type == "truth":
            pipeline_fac = pipeline_factory.Thruth

        elif container.type == "comp":
            pipeline_fac = pipeline_factory.Comp

        for target_pipe in container.pipe:
            log_print(f"\n\n[CREATING PIPELINE] {target_pipe}")
            pipe = pipeline_fac.build_pipeline(container, target_pipe)

            log_print(f"\n[EXECUTING]")
            success, times = pipe.execute()
            if not success:
                log_print(f"[FAILURE]")
                log_print(f"--- TEST FAILED AT EXECUTION---")
                container.success[target_pipe] = 1
                continue

            log_print(f"\n[COMPARING RESULTS]")
            success = pipe.compare()
            if not success:
                log_print(f"[FAILURE]")
                log_print(f"--- TEST FAILED AT COMPARISON ---")
                container.success[target_pipe] = 2
            container.times[target_pipe] = times

        log_print("\n[CONTAINER SUMMARY]")
        run_map = ["success", "crash", "failed"]
        for item in container.success.items():
            pipe, success = item
            log_print(f"\t{pipe}:\t{run_map[success]}")
        
        if container.type == "comp":
            log_print(f"\n[PERFORMANCE SUMMARY]")
            for item in container.times.items():
                pipe, times = item 
                if container.success[pipe] == 0:
                    log_print(f"\t{pipe}:{time.compare_and_format(*times)}")
        console.export_sublog(container)

    run_map = ["success", "crash", "failed"]
    log_print("\n--- [CI SUMMARY] ---")
    for container in containers:
        log_print(f"\t[{container.name}]")
        for item in container.success.items():
            pipe, success = item
            log_print(f"\t\t{pipe}:\t{run_map[success]}")

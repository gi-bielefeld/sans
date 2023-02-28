"""
This file implements the parsing of run
scripts
"""

import json

from . import path
from . import container
from .console import log_print

__VALID_RUN_TYPES = ["truth", "comp", "perf"]
__VALID_DATA_TYPES = ["synthetic", "comp", "external"]

__VALID_PIPES = {"truth": ["input", "nrev", "amino"], "comp": ["input", "nrev", "amino", "trans", "iupac", "window"]}


def __throw_parsing_error(file, info):
    """
    This method throws an error, logs it and exits the application
    :param file: The file that was checked, when the error occurred
    :param info: Information to print
    :return:
    """
    log_print(f"[SCRIPT_PARSER] ERROR: Failed parsing script file {file}\n\t{info}")
    exit(1)


def __validate_data(file: str, run_name: str, container_type: str, data: dict):
    """
    This method checks whether the given definition of the data to use is valid
    :param container_type: The container type
    :param data: The data definition
    :return: None, Throws parsing error if invalid
    """
    if container_type == "truth" and ("source" in data and data["source"] != "synthetic"):
        __throw_parsing_error(file, f"Run {run_name}: 'truth' type testing requires synthetic data")

    if data["source"] in ["synthetic", "comp"]:
        # data.seq_num
        if "seq_num" not in data:
            __throw_parsing_error(file, f"Run {run_name}:  Missing definition of data.seq_num")
        elif type(data["seq_num"]) != int:
            __throw_parsing_error(file, f"Run {run_name}:  Invalid definition of data.seq_num: {data['seq_num']}")
        elif data["seq_num"] < 2:
            __throw_parsing_error(file, f"Run {run_name}:  Invalid definition of data.seq_num: {data['seq_num']} < 2")

        # data.seq_len
        elif "seq_len" not in data:
            __throw_parsing_error(file, f"Run {run_name}:  Missing definition of data.seq_len")
        elif type(data["seq_len"]) != int:
            __throw_parsing_error(file, f"Run {run_name}:  Invalid definition of data.seq_len: {data['seq_len']}")

        # data.m_rate
        elif "m_rate" not in data:
            __throw_parsing_error(file, f"Run {run_name}:  Missing definition of data.m_rate: {data['m_rate']}")
        elif type(data["m_rate"]) != int:
            __throw_parsing_error(file, f"Run {run_name}:  Invalid definition of data.mrate: {data['m_rate']}")
        elif 0 > data["m_rate"] or data["m_rate"] > 10:
            __throw_parsing_error(file, f"Run {run_name}:  Invalid definition of data.mrate: 0 < {data['m_rate']} < 10")

    elif data["source"] == "external":
        if "path" not in data:
            __throw_parsing_error(file, f"Run {run_name}:  Missing definition of data.path")
        elif not path.exists(data["path"]) or path.is_dir(data[path]):
            __throw_parsing_error(file, f"Run {run_name}:  External data file not found {data['path']}")

    elif "source" in data:
        __throw_parsing_error(file, f"Run {run_name}: Invalid definition of data source: {data['source']}")


def __validate_sans(file: str, run_name: str, container_type: str, sans: dict):
    # sans.k
    if "k" not in sans:
        __throw_parsing_error(file, f"Run {run_name}, Missing definition of sans.k")
    elif type(sans["k"]) != int:
        __throw_parsing_error(file, f"Run {run_name}, Invalid definition of sans.k: {sans['k']}")

    # TRUTH CONTAINERS
    if container_type == "truth":
        if "target_a" not in sans:
            __throw_parsing_error(file, f"Run {run_name}, Missing definition of sans.target_a")
        elif type(sans["target_a"]) != str:
            __throw_parsing_error(file, f"Run {run_name}, Invalid definition of sans.target_a: {sans['target_a']}")
        if "target_b" in sans:
            __throw_parsing_error(file, f"Run {run_name}, Truth test does not support two targets sans.target_b: {sans['target_b']}")
    # COMP CONTAINERS
    elif container_type == "comp":
        if "target_a" not in sans:
            __throw_parsing_error(file, f"Run {run_name}, Missing definition of sans.target_a")
        elif type(sans["target_a"]) != str:
            __throw_parsing_error(file, f"Run {run_name}, Invalid definition of sans.target_a: {sans['target_a']}")

        if "target_b" not in sans:
            __throw_parsing_error(file, f"Run {run_name}, Missing definition of sans.target_b")
        elif type(sans["target_b"]) != str:
            __throw_parsing_error(file, f"Run {run_name}, Invalid definition of sans.target_b: {sans['target_b']}")



def __validate_pipe(file: str, run_name: str, container_type: str, pipe: list):
    """
    This method checks if the given test pipes are supported by the container type
    :param file: The current json file
    :param run_name: The name of the checked run
    :param container_type: The target container type
    :param pipe: The target pipelines
    :return: None, Throws parsing error
    """
    for element in pipe:
        if element not in __VALID_PIPES[container_type]:
            __throw_parsing_error(file, f"Run {run_name}, Invalid pipe '{element}' for container {container_type}")


def parse_run_script(run_script):
    log_print(f"[SCRIPT_PARSER] Parsing run script: {run_script}")
    # Check if the given file exists:
    if not path.exists(run_script) or path.is_dir(run_script):
        __throw_parsing_error(run_script, f"Run file not found, at: {run_script}")

    # Read lines from the file
    file = open(run_script)
    content = file.read()
    file.close()

    # Parse
    runs = json.loads(content)

    # Container aggregator
    containers = []

    for run_name in runs.keys():
        elements = runs[run_name]

        # Check the run type
        check_type = True
        if "type" not in elements:
            __throw_parsing_error(file, f"Run {run_name}: Missing 'type' definition")

        elif elements["type"] not in __VALID_RUN_TYPES:
            __throw_parsing_error(file, f"Run {run_name}: Invalid 'type' definition: {elements['type']}")

        if "data" not in elements:
            __throw_parsing_error(file, f"Run {run_name}: Missing 'data' definition")

        if "sans" not in elements:
            __throw_parsing_error(file, f"Run {run_name}: Missing 'sans' definition")

        if "pipe" not in elements:
            __throw_parsing_error(file, f"Run {run_name}: Missing 'pipe' definition")

        run_type = elements["type"]
        data = elements["data"]
        sans = elements["sans"]
        pipe = elements["pipe"]

        __validate_data(run_script, run_name, run_type, data)
        __validate_sans(run_script, run_name, run_type, sans)
        __validate_pipe(run_script, run_name, run_type, pipe)
        # Build containers based on the run_type
        containers.append(container.Container(name=run_name, **elements))
    return containers

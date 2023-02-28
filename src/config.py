
from . import path
from . console import log_print

import json


def throw_config_error(msg):
    log_print(f"[CONFIG] ERROR:  \n{msg}")
    exit(1)


class Config:

    using_default_root = None
    root_path = None
    seq_path = None
    data_path = None
    log_path = None

    @classmethod
    def init(cls):
        """
        This method loads the global configurations from the config file
        :return: None
        """
        log_print("[CONFIG] Loading Config")
        config_path = path.nodes_to_dir_path([path.ROOT, "config.json"])
        if not path.exists(config_path):
            throw_config_error("Config file not found")

        config_file = open(config_path)
        content = config_file.read()
        config_file.close()

        config = json.loads(content)
        if "ROOT" not in config:
            throw_config_error("Missing definition of ROOT")

        root = config["ROOT"]
        if root != "DEFAULT_ROOT":
            cls.using_default_root = False
            cls.root_path = root
        else:
            cls.root_path = path.nodes_to_dir_path([path.ROOT, root])

        if "SEQ_PATH" not in config:
            throw_config_error("Missing definition of SEQ_PATH")
        seq_path = config["SEQ_PATH"]
        cls.seq_path = path.nodes_to_dir_path([cls.root_path, seq_path])

        if "DATA_PATH" not in config:
            throw_config_error("Missing definition of DATA_PATH")
        data_path = config["DATA_PATH"]
        cls.data_path = path.nodes_to_dir_path([cls.root_path, data_path])

        if "LOG_PATH" not in config:
            throw_config_error("Missing definition of LOG_PATH")

        log_path = config["LOG_PATH"]
        cls.log_path = path.nodes_to_dir_path([cls.root_path, log_path])

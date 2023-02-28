"""
This file implements basic path
handling and global paths
"""

import os
from pathlib import Path


ROOT = "/".join(str(Path(__file__).parent.parent.absolute()).split("\\"))


def __clean_path_nodes(path_nodes: list):
    """
    This method ensures none of the nodes
    starts or ends with a separator
    :param path_nodes:
    :return: list: clean_nodes
    """
    clean_nodes = []
    for node in path_nodes:
        clean_node = node
        # remove the last character if it is a seqarator
        if clean_node.endswith("/") or clean_node.endswith("\\"):
            clean_node = clean_node[:-1]
        # add the node if it is not empty
        clean_nodes.append(clean_node)
    return clean_nodes


def nodes_to_file_path(nodes):
    """
    This method merges a set of directories into a path
    :param nodes: The nodes to merge
    :return: str: file_path
    """
    return "/".join(__clean_path_nodes(nodes))


def nodes_to_dir_path(nodes):
    """
    This method merges a set of directories into
    a path in directory format
    :param nodes: The nodes to merge
    :return: str: dir_path
    """
    return "/".join(__clean_path_nodes(nodes))


def exists(path):
    """
    This method checks whether the given path exists in the system
    :param path: The path to check
    :return: bool: Exists
    """
    return os.path.exists(path)


def is_dir(path):
    """
    This method checks whether the given path points to an existing directory
    :param path: The path to check
    :return: bool: True only if the path exists and points to a directory
    """
    if not exists(path):
        return False
    return os.path.isdir(path)


def touch_file(path):
    """
    This method touches a file in the file system
    :return: None
    """
    # Catch file already exists
    if exists(path):
        print(f"[Error: ] Trying to create existing path {path}")
        exit(1)

    # Catch path is in directoy format
    elif is_dir(path):
        print(f"[Error: ] Trying to create a directory as file {path}")
        exit(1)

    # Create the file
    else:
        f = open(path, "w+")
        f.close()


def touch_directory(path):
    # Catch file already exists
    if exists(path) and is_dir(path):
        print(f"[Error: ] Trying to create existing path {path}")
        exit(1)

    # Create the file
    else:
        os.makedirs(path)

"""
This file implement the container for
the runs meta information.
@Author Fabian Kolesch
"""


class Data:
    """
    This class stores information about the data to use for a run
    """
    def __init__(self, source: str, path=None, seq_num: int = 0, seq_len: int = 0, m_rate: float = 0):
        self.source = source
        self.external_path = path
        self.seq_num = seq_num
        self.seq_len = seq_len
        self.m_rate = m_rate


class Sans:
    """
    This class stores sans parameters
    """
    def __init__(self, k, target_a, target_b=None):
        self.k = k
        self.target_a = target_a
        self.target_b = target_b


# Template for run container
class Container:
    """
    This class stores information about the run that is to be executed
    """
    def __init__(self, name: str, type: str, data: dict, sans: dict, pipe: list):
        self.name = name
        self.type = type
        self.data = Data(**data)
        self.sans = Sans(**sans)
        self.pipe = pipe

        self.seq_dir = None
        self.data_dir = None
        self.log_dir = None

        self.success = dict()
        for element in self.pipe:
            self.success[element] = 0

    def get_info(self):
        info = ""
        info += f"name:\t\t{self.name}\ntype:\t\t{self.type}\n"
        info += f"data:\t\t{vars(self.data)}\nsans:\t\t{vars(self.sans)}\npipe:\t\t{self.pipe}\n"
        return info
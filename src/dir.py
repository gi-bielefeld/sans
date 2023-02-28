
from .config import Config
from .console import log_print
from .container import Container
from . import path


class Seq:
    def __init__(self, container: Container):
        log_print(f"[SEQ] Touching")
        seq_root = Config.seq_path
        if not path.exists(seq_root):
            path.touch_directory(seq_root)
        if container.type == "truth":
            set_header = "_".join(['Tru',
                                   'k' + str(container.sans.k),
                                   'n' + str(container.data.seq_num),
                                   'l' + str(container.data.seq_len),
                                   'm' + str(container.data.m_rate),
                                   ])
        else:
            set_header = "_".join(['Com',
                                   'k' + str(container.sans.k),
                                   'n' + str(container.data.seq_num),
                                   'l' + str(container.data.seq_len),
                                   'm' + str(container.data.m_rate),
                                   ])

        set_path = path.nodes_to_dir_path([seq_root, set_header])
        if not path.exists(set_path) or not path.is_dir(set_path):
            path.touch_directory(set_path)

        fa_path = path.nodes_to_dir_path([set_path, "fa"])
        if not path.exists(fa_path) or not path.is_dir(fa_path):
            path.touch_directory(fa_path)
        log_print(f"\t{set_path}")

        self.set_path = set_path
        self.fa_path = fa_path

        # TRUTH fields
        # Input
        self.dna_list_path = None
        self.dna_splits_path = None
        # NoReverse
        self.dna_nrev_list_path = None
        self.dna_nrev_splits_path = None
        # Amino
        self.amino_list_path = None
        self.amino_splits_path = None

        # Set this to be the containers seq_dir
        container.seq_dir = self


def _throw_data_error(msg):
    log_print(f"[DATA] ERROR: {msg}")
    exit(1)


class RunData:
    def __init__(self, container):
        log_print(f"[RUN_DATA] Touching")
        data_path = Config.data_path
        if not path.exists(data_path):
            path.touch_directory(data_path)

        container_path = path.nodes_to_dir_path([data_path, container.name])
        if not path.exists(container_path) or not path.is_dir(container_path):
            path.touch_directory(container_path)

        log_print(f"\t{container_path}")

        self.container_path = container_path

        # TRUTH COMPATIBLE PIPELINES
        # INPUT
        self.input_target_a_splits_path = path.nodes_to_file_path([self.container_path, "INPUT_target_a_splits.txt"])
        self.input_target_b_splits_path = path.nodes_to_file_path([self.container_path, "INPUT_target_b_splits.txt"])
        self.input_diff_path = path.nodes_to_file_path([self.container_path, "INPUT_splits_diff.txt"])

        # NREV
        self.nrev_target_a_splits_path = path.nodes_to_file_path([self.container_path, "NREV_target_a_splits.txt"])
        self.nrev_target_b_splits_path = path.nodes_to_file_path([self.container_path, "NREV_target_b_splits.txt"])
        self.nrev_diff_path = path.nodes_to_file_path([self.container_path, "NREV_splits_diff.txt"])

        # AMINO
        self.amino_target_a_splits_path = path.nodes_to_file_path([self.container_path, "AMINO_target_a_splits.txt"])
        self.amino_target_b_splits_path = path.nodes_to_file_path([self.container_path, "AMINO_target_b_splits.txt"])
        self.amino_diff_path = path.nodes_to_file_path([self.container_path, "AMINO_splits_diff.txt"])
        
        # OTHER PIPELINES
        # TRANS
        self.trans_target_a_splits_path = path.nodes_to_file_path([self.container_path, "TRANS_target_a_splits.txt"])
        self.trans_target_b_splits_path = path.nodes_to_file_path([self.container_path, "TRANS_target_b_splits.txt"])
        self.trans_diff_path = path.nodes_to_file_path([self.container_path, "TRANS_splits_diff.txt"])
        # IUPAC
        self.iupac_target_a_splits_path = path.nodes_to_file_path([self.container_path, "IUPAC_target_a_splits.txt"])
        self.iupac_target_b_splits_path = path.nodes_to_file_path([self.container_path, "IUPAC_target_b_splits.txt"])
        self.iupac_diff_path = path.nodes_to_file_path([self.container_path, "IUPAC_splits_diff.txt"])       
        # WINDOW
        self.window_target_a_splits_path = path.nodes_to_file_path([self.container_path, "WINDOW_target_a_splits.txt"])
        self.window_target_b_splits_path = path.nodes_to_file_path([self.container_path, "WINDOW_target_b_splits.txt"])
        self.window_diff_path = path.nodes_to_file_path([self.container_path, "WINDOW_splits_diff.txt"])   

        container.data_dir = self


def _throw_test_log_error(msg):
    log_print(f"[TEST_LOG] ERROR: {msg}")
    exit(1)


class TestLog:

    def __init__(self, container):
        log_print(f"[TEST_LOG] Touching")
        test_log_path = Config.log_path
        if not path.exists(test_log_path):
            path.touch_directory(test_log_path)

        container_path = path.nodes_to_dir_path([test_log_path, container.name])
        if not path.exists(container_path) or not path.is_dir(container_path):
            path.touch_directory(container_path)
        container_path = container_path

        log_print(f"\t{container_path}")

        self.container_path = container_path
        self.log_path = path.nodes_to_file_path([container_path, "log.txt"])

        container.log_dir = self

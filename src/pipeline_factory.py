

import subprocess

from .container import Container
from .console import log_print


class _Executor:
    def __init__(self, commands):
        self.commands = commands

    def run(self):
        success = True
        for command in self.commands:
            log_print(f"[EXECUTING COMMAND]:\n{command}")
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            retval = proc.wait()
            success = retval == 0
            out = out.decode("UTF-8")
            err = err.decode("UTF-8")
            log_print(f"[STDOUT] -------")
            log_print(out)
            log_print(f"[CERR] (fatal={not success})-------")
            log_print(err)
            log_print(f"... DONE")
            if not success:
                break
        return success


class _Comparator:
    def __init__(self, ref_list, path_a, path_b):
        self.ref_list = ref_list
        self.path_a = path_a
        self.path_b = path_b

    def run(self):
        log_print(f"[COMPARING]")
        log_print(f"\tthis one: {self.path_a}")
        log_print(f"\tthat one: {self.path_b}")
        # --- [Load list]
        list_file = open(self.ref_list, 'r')
        lines = list_file.readlines()
        list_file.close()
        seq_names = [line.strip("\n") for line in lines]
        seq_name_map = dict()
        for i in range(len(seq_names)):
            seq_name_map[seq_names[i]] = i

        # --- [Read splits] ---
        # First file
        file_a = open(self.path_a, 'r')
        lines_a = file_a.readlines()
        file_a.close()
        splits_a = dict()
        for split_line in lines_a:
            split_line = split_line.strip("\n")
            fields = split_line.split("\t")
            # The split weight
            weight = float(fields[0])
            # The file names
            split_seq_names = fields[1:]
            # Translate file names -> int(sortable) -> string(hashable)
            colors = [seq_name_map[split_seq_name] for split_seq_name in split_seq_names]
            colors = list(sorted(colors))
            color_string = ','.join([str(color) for color in colors])
            # Hash the split
            if color_string not in splits_a:
                splits_a[color_string] = set()
            splits_a[color_string].add(weight)

        # Other file
        file_b = open(self.path_b, 'r')
        lines_b = file_b.readlines()
        file_b.close()
        splits_b = dict()
        for split_line in lines_b:
            split_line = split_line.strip("\n")
            fields = split_line.split("\t")
            # The split weight
            weight = float(fields[0])
            # The file names
            split_seq_names = fields[1:]
            # Translate file names -> int(sortable) -> string(hashable)
            colors = [seq_name_map[split_seq_name] for split_seq_name in split_seq_names]
            colors = list(sorted(colors))
            color_string = ','.join([str(color) for color in colors])
            # Hash the split
            if color_string not in splits_b:
                splits_b[color_string] = set()
            splits_b[color_string].add(weight)

        for item in splits_a.items():
            colors, weight = item
            if colors not in splits_b:
                return False
            # elif weight not in splits_b[colors]:
            #    return False

        for item in splits_b.items():
            colors, weight = item
            if colors not in splits_a:
                return False
            # elif weight not in splits_b[colors]:
            #    return False
        return True


class _Pipeline:
    def __init__(self, executor, comparator):
        self._executor = executor
        self._comparator = comparator

    def execute(self):
        return self._executor.run()

    def compare(self):
        return self._comparator.run()


class Thruth:
    @classmethod
    def _input(cls, container: Container):
        input_command = [container.sans.target_a,
                         f" -i {container.seq_dir.dna_list_path}",
                         f" -o {container.data_dir.input_target_a_splits_path}"
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command = "".join(input_command)

        ref_list = container.seq_dir.dna_list_path
        ref_splits = container.seq_dir.dna_splits_path
        comp_splits = container.data_dir.input_target_a_splits_path

        executor = _Executor([input_command])
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def _nrev(cls, container: Container):
        input_command = [container.sans.target_a,
                         f" -i {container.seq_dir.dna_nrev_list_path}",
                         f" -o {container.data_dir.nrev_target_a_splits_path}"
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command = "".join(input_command)

        ref_list = container.seq_dir.dna_list_path
        ref_splits = container.seq_dir.dna_splits_path
        comp_splits = container.data_dir.nrev_target_a_splits_path

        executor = _Executor([input_command])
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def _amino(cls, container: Container):
        input_command = [container.sans.target_a,
                         f" -i {container.seq_dir.amino_list_path}",
                         f" -o {container.data_dir.amino_target_a_splits_path}",
                         f" -a",
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command = "".join(input_command)

        ref_list = container.seq_dir.amino_list_path
        ref_splits = container.seq_dir.amino_splits_path
        comp_splits = container.data_dir.amino_target_a_splits_path

        executor = _Executor([input_command])
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def build_pipeline(cls, container: Container, target_pipe):
        pipe_map = {"input": cls._input, "nrev": cls._nrev, "amino": cls._amino}
        builder = pipe_map[target_pipe]
        pipeline = builder(container)
        return pipeline


class Comp:
    @classmethod
    def _input(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                         f" -i {container.seq_dir.dna_list_path}",
                         f" -o {container.data_dir.input_target_a_splits_path}"
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)
        # Run command for target_b
        input_command_b = [container.sans.target_b,
                         f" -i {container.seq_dir.dna_list_path}",
                         f" -o {container.data_dir.input_target_b_splits_path}"
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_list = container.seq_dir.dna_list_path
        ref_splits = container.data_dir.input_target_a_splits_path
        comp_splits = container.data_dir.input_target_b_splits_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def _nrev(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.dna_nrev_list_path}",
                           f" -o {container.data_dir.nrev_target_a_splits_path}",
                           f" -m geom2",
                           f" -n"]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)
        # Run command for target_b
        input_command_b = [container.sans.target_b,
                           f" -i {container.seq_dir.dna_nrev_list_path}",
                           f" -o {container.data_dir.nrev_target_b_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -n"]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_list = container.seq_dir.dna_nrev_list_path
        ref_splits = container.data_dir.nrev_target_a_splits_path
        comp_splits = container.data_dir.nrev_target_b_splits_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def _amino(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.amino_list_path}",
                           f" -o {container.data_dir.amino_target_a_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -a"]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)

        # Run command for target_b
        input_command_b = [container.sans.target_b,
                           f" -i {container.seq_dir.amino_list_path}",
                           f" -o {container.data_dir.amino_target_b_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -a"]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_list = container.seq_dir.amino_list_path
        ref_splits = container.data_dir.amino_target_a_splits_path
        comp_splits = container.data_dir.amino_target_b_splits_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_list, ref_splits, comp_splits)
        return _Pipeline(executor, comparator)

    @classmethod
    def build_pipeline(cls, container: Container, target_pipe):
        pipe_map = {"input": cls._input, "nrev": cls._nrev, "amino": cls._amino}
        builder = pipe_map[target_pipe]
        pipeline = builder(container)
        return pipeline

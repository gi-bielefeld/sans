

import subprocess

from .container import Container
from .console import log_print

from .time import Timer


class _Executor:
    def __init__(self, commands):
        self.commands = commands

    def run(self):
        success = True
        timer = Timer()
        times = []
        for command in self.commands:
            log_print(f"[EXECUTING COMMAND]:\n{command}")
            timer.start()
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            retval = proc.wait()
            timer.stop()
            success = retval == 0
            times.append(timer.delta)
            out = out.decode("UTF-8")
            err = err.decode("UTF-8")
            log_print(f"[STDOUT] -------")
            log_print(out)
            log_print(f"[CERR] (fatal={not success})-------")
            log_print(err)
            log_print(f"... DONE")
            if not success:
                break
        return success, times


class _Comparator:
    def __init__(self, path_a, path_b, path_diff):
        self.path_a = path_a
        self.path_b = path_b
        self.diff_path = path_diff

    def run(self):
        commands = []
        commands.append(f"sort -u {self.path_a} >> {self.path_a}.sorted")
        commands.append(f"sort -u {self.path_b} >> {self.path_b}.sorted")
        commands.append(f"comm -3 {self.path_a}.sorted {self.path_b}.sorted  >> {self.diff_path}")
        for command in commands:
            log_print(f"[EXECUTING COMMAND]:\n{command}")
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            retval = proc.wait()
        
        diff_file = open(self.diff_path, 'r')
        lines = diff_file.readlines()
        diff_file.close()
        return not lines


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

        ref_splits = container.seq_dir.dna_splits_path
        comp_splits = container.data_dir.input_target_a_splits_path
        diff_path = container.data_dir.input_diff_path

        executor = _Executor([input_command])
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
        return _Pipeline(executor, comparator)

    @classmethod
    def _nrev(cls, container: Container):
        input_command = [container.sans.target_a,
                         f" -i {container.seq_dir.dna_nrev_list_path}",
                         f" -o {container.data_dir.nrev_target_a_splits_path}",
                         f" -k {container.sans.k}",
                         f" -m geom2"]
        input_command = "".join(input_command)

        ref_splits = container.seq_dir.dna_nrev_splits_path
        comp_splits = container.data_dir.nrev_target_a_splits_path
        diff_path = container.data_dir.nrev_diff_path 

        executor = _Executor([input_command])
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
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

        ref_splits = container.seq_dir.amino_splits_path
        comp_splits = container.data_dir.amino_target_a_splits_path
        diff_path = container.data_dir.amino_diff_path

        executor = _Executor([input_command])
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
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

        ref_splits = container.data_dir.input_target_a_splits_path
        comp_splits = container.data_dir.input_target_b_splits_path
        diff = container.data_dir.input_diff_path        
        
        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff)
        return _Pipeline(executor, comparator)

    @classmethod
    def _nrev(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.dna_nrev_list_path}",
                           f" -o {container.data_dir.nrev_target_a_splits_path}",
                           f" -k {container.sans.k}",
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

        ref_splits = container.data_dir.nrev_target_a_splits_path
        comp_splits = container.data_dir.nrev_target_b_splits_path
        diff_path = container.data_dir.nrev_diff_path           
    
        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
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

        ref_splits = container.data_dir.amino_target_a_splits_path
        comp_splits = container.data_dir.amino_target_b_splits_path
        diff_path = container.data_dir.amino_diff_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
        return _Pipeline(executor, comparator)

    @classmethod
    def _trans(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.trans_target_a_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -c"
                          ]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)

        # Run command for target_b
        input_command_b = [container.sans.target_b,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.trans_target_b_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -c"
                          ]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_splits = container.data_dir.trans_target_a_splits_path
        comp_splits = container.data_dir.trans_target_b_splits_path
        diff_path = container.data_dir.trans_diff_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
        return _Pipeline(executor, comparator)

    @classmethod
    def _iupac(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.iupac_target_a_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -x 11"
                          ]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)

        # Run command for target_b
        input_command_b = [container.sans.target_b,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.iupac_target_b_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -x 11"
                          ]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_splits = container.data_dir.iupac_target_a_splits_path
        comp_splits = container.data_dir.iupac_target_b_splits_path
        diff_path = container.data_dir.iupac_diff_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
        return _Pipeline(executor, comparator)

    @classmethod
    def _window(cls, container: Container):
        commands = []
        # Run command for target_a
        input_command_a = [container.sans.target_a,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.window_target_a_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -w 11"
                          ]
        input_command_a = "".join(input_command_a)
        commands.append(input_command_a)

        # Run command for target_b
        input_command_b = [container.sans.target_b,
                           f" -i {container.seq_dir.dna_list_path}",
                           f" -o {container.data_dir.window_target_b_splits_path}",
                           f" -k {container.sans.k}",
                           f" -m geom2",
                           f" -w 11"
                          ]
        input_command_b = "".join(input_command_b)
        commands.append(input_command_b)

        ref_splits = container.data_dir.window_target_a_splits_path
        comp_splits = container.data_dir.window_target_b_splits_path
        diff_path = container.data_dir.window_diff_path

        executor = _Executor(commands)
        comparator = _Comparator(ref_splits, comp_splits, diff_path)
        return _Pipeline(executor, comparator)


    @classmethod
    def build_pipeline(cls, container: Container, target_pipe):
        pipe_map = {"input": cls._input, "nrev": cls._nrev, "amino": cls._amino, "trans": cls._trans, "iupac": cls._iupac, "window": cls._window}
        builder = pipe_map[target_pipe]
        pipeline = builder(container)
        return pipeline

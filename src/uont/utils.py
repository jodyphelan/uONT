import logging
import subprocess as sp
import re
from shutil import which
import os
import tempfile
import shutil
import functools
import inspect
import time
from .types import FullPath

def run_cmd(cmd: str, desc=None, log: str=None, exit_on_error: bool=True) -> sp.CompletedProcess:
    if desc:
        logging.info(desc)
    processed_cmd = cmd.replace("&&","XX")
    programs = set([x.strip().split()[0] for x in re.split("[|&;]",processed_cmd.strip()) if x!=""])
    missing = [p for p in programs if which(p)==False]
    if len(missing)>0:
        raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
    logging.debug(f"Running command: {cmd}")
    cmd = "/bin/bash -c set -o pipefail; " + cmd
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,stderr=output,stdout=output)
    if result.returncode != 0:
        logging.error(result.stderr.decode("utf-8"))
        if exit_on_error:
            raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,result.stderr.decode()))
    return result

def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.info(f"Execution time for {func.__name__}: {elapsed_time:.2f} seconds")
        return result
    return wrapper

def run_in_tempdir(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        try:
            # filter out any arguments in kwargs that are already in args to avoid duplication
            sig = inspect.signature(func)

            required_arg_names = [
                param.name for param in sig.parameters.values()
                if param.default == inspect.Parameter.empty
                and param.name != 'kwargs'
            ]

            args_dict = {}
            for i, param_name in enumerate(required_arg_names):
                if i < len(args):
                    args_dict[param_name] = args[i]
                    if param_name in kwargs:
                        del kwargs[param_name]
                elif param_name in kwargs:
                    args_dict[param_name] = kwargs[param_name]
                    del kwargs[param_name]


            args = tuple(args_dict.values())
            
            # convert any FullPath arguments to absolute paths
            sig = inspect.signature(func)
            for param in sig.parameters.values():
                arg_type = param.annotation if param.annotation is not param.empty else str

                if arg_type == FullPath:
                    arg_index = list(sig.parameters).index(param.name)

                    if arg_index < len(args):
                        args = list(args)
                        args[arg_index] = os.path.abspath(args[arg_index])
                        args = tuple(args)

                    elif param.name in kwargs:
                        if kwargs[param.name] is not None:
                            kwargs[param.name] = os.path.abspath(kwargs[param.name])
                    else:
                        logging.debug(f"Argument {param.name} not found in args or kwargs, skipping path conversion")

            logging.debug(f"Arguments for {func.__name__}: args={args}, kwargs={kwargs}")
            logging.debug(f"Running {func.__name__} in temporary directory: {tmpdir}")
            os.chdir(tmpdir)

            return func(*args, tmp_dir=tmpdir, **kwargs)
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir, ignore_errors=True)
    return wrapper

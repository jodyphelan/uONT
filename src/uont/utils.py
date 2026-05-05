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


DEFAULT_CLI_DEPENDENCIES = [
    "autocycler",
    "bcftools",
    "chopper",
    "dnaapler",
    "dorado",
    "filtlong",
    "flye",
    "hostile",
    "lrge",
    "medaka",
    "miniasm",
    "minimap2",
    "pigz",
    "porechop_abi",
    "rmlst",
    "samtools",
    "seqkit",
]

def get_software_version(tool: str) -> str:
    cmds = {
        "autocycler": "autocycler --version",
        "bcftools": "bcftools --version",
        "chopper": "chopper --version",
        "dnaapler": "dnaapler --version",
        "dorado": "dorado --version",
        "filtlong": "filtlong --version",
        "flye": "flye --version",
        "hostile": "hostile --version",
        "lrge": "lrge --version",
        "medaka": "medaka --version",
        "miniasm": "miniasm -V",
        "minimap2": "minimap2 --version",
        "pigz": "pigz --version",
        "porechop_abi": "porechop_abi --version",
        "rmlst": "rmlst --version",
        "samtools": "samtools --version",
        "seqkit": "seqkit version",
    }
    regex = {
        "autocycler": r"autocycler\s+v?([0-9][0-9A-Za-z_.-]*)",
        "bcftools": r"bcftools\s+([0-9][0-9A-Za-z_.-]*)",
        "chopper": r"chopper\s+v?([0-9][0-9A-Za-z_.-]*)",
        "dnaapler": r"dnaapler\s+v?([0-9][0-9A-Za-z_.-]*)",
        "dorado": r"dorado\s+v?([0-9][0-9A-Za-z_.-]*)",
        "filtlong": r"filtlong\s+v?([0-9][0-9A-Za-z_.-]*)",
        "flye": r"flye\s+v?([0-9][0-9A-Za-z_.-]*)",
        "hostile": r"hostile\s+v?([0-9][0-9A-Za-z_.-]*)",
        "lrge": r"lrge\s+v?([0-9][0-9A-Za-z_.-]*)",
        "medaka": r"medaka\s+v?([0-9][0-9A-Za-z_.-]*)",
        "miniasm": r"miniasm\s+v?([0-9][0-9A-Za-z_.-]*)",
        "minimap2": r"([0-9][0-9A-Za-z_.-]*)",
        "pigz": r"pigz\s+([0-9][0-9A-Za-z_.-]*)",
        "porechop_abi": r"porechop(?:_abi)?\s+v?([0-9][0-9A-Za-z_.-]*)",
        "rmlst": r"rmlst\s+v?([0-9][0-9A-Za-z_.-]*)",
        "samtools": r"samtools\s+([0-9][0-9A-Za-z_.-]*)",
        "seqkit": r"seqkit\s+v?([0-9][0-9A-Za-z_.-]*)",
    }

    if tool not in cmds:
        raise ValueError(f"No version command configured for tool: {tool}")

    x = sp.run(cmds[tool], shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    text = x.stdout.decode("utf-8", errors="ignore") + x.stderr.decode("utf-8", errors="ignore")
    match = re.search(regex[tool], text, flags=re.IGNORECASE)
    if not match:
        # Fallback: try to extract the first semantic-looking version from output.
        fallback = re.search(r"([0-9]+\.[0-9]+(?:\.[0-9A-Za-z_.-]+)?)", text)
        if fallback:
            return fallback.group(1)
        raise ValueError(f"Could not parse version for tool '{tool}' from output:\n{text}")
    return match.group(1)

def check_cli_dependencies(programs: list[str] | None = None) -> tuple[list[str], list[str]]:
    """Return available and missing CLI dependencies.

    Args:
        programs (list[str] | None): CLI programs to check. If None, checks
            ``DEFAULT_CLI_DEPENDENCIES``.

    Returns:
        tuple[list[str], list[str]]: (available, missing)
    """
    programs_to_check = programs if programs is not None else DEFAULT_CLI_DEPENDENCIES

    # Preserve order while avoiding duplicate checks.
    seen = set()
    unique_programs = []
    for program in programs_to_check:
        if program not in seen:
            seen.add(program)
            unique_programs.append(program)

    available: list[str] = []
    missing: list[str] = []

    for program in unique_programs:
        if which(program) is None:
            missing.append(program)
            print(f"🔴 {program}: not found")
            continue

        available.append(program)
        try:
            version = get_software_version(program)
            print(f"🟢 {program}: {version}")
        except Exception:
            print(f"🟠 {program}: found, version unknown")

    return available, missing

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

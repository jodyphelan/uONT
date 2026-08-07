import logging
from enum import Enum
import subprocess as sp
import re
import signal
from shutil import which
import os
import tempfile
import shutil
import functools
import inspect
import time
from typing import Optional
from .types import FullPath
from dataclasses import fields
import uont

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
    "NanoPlot",
    "plassembler",
    "myloasm",
    "nextDenovo",
    "nextPolish",
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
        "NanoPlot": "NanoPlot --version",
        "plassembler": "plassembler --version",
        "myloasm": "myloasm --version",
        "nextDenovo": "nextDenovo --version",
        "nextPolish": "nextPolish --version",
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
        "NanoPlot": r"NanoPlot\s+v?([0-9][0-9A-Za-z_.-]*)",
        "plassembler": r"plassembler\s+v?([0-9][0-9A-Za-z_.-]*)",
        "myloasm": r"myloasm\s+v?([0-9][0-9A-Za-z_.-]*)",
        "nextDenovo": r"nextDenovo\s+v?([0-9][0-9A-Za-z_.-]*)",
        "nextPolish": r"nextPolish\s+v?([0-9][0-9A-Za-z_.-]*)",
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

def run_cmd(
    cmd: str,
    desc=None,
    log: Optional[str] = None,
    exit_on_error: bool=True,
    timeout: Optional[float] = None,
) -> sp.CompletedProcess:
    if desc:
        logging.info(desc)
    processed_cmd = cmd.replace("&&","XX")
    programs = set([x.strip().split()[0] for x in re.split("[|&;]",processed_cmd.strip()) if x!=""])
    missing = [p for p in programs if which(p)==False]
    if len(missing)>0:
        raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
    logging.debug(f"Running command: {cmd}")
    log_handle = open(log,"w") if log else None
    output = log_handle if log_handle else sp.PIPE
    process = None
    try:
        process = sp.Popen(
            ["/bin/bash", "-o", "pipefail", "-c", cmd],
            stdout=output,
            stderr=output,
            start_new_session=True,
        )
        stdout, stderr = process.communicate(timeout=timeout)
        result = sp.CompletedProcess(process.args, process.returncode, stdout, stderr)
    except sp.TimeoutExpired as exc:
        if process is not None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                stdout, stderr = process.communicate(timeout=5)
            except sp.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                stdout, stderr = process.communicate()
        else:
            stdout = exc.stdout
            stderr = exc.stderr
        timeout_message = f"Command timed out after {timeout} seconds:\n{cmd}"
        logging.error(timeout_message)
        if exit_on_error:
            raise TimeoutError(timeout_message) from exc
        result = sp.CompletedProcess(cmd, 124, stdout, stderr)
    finally:
        if log_handle is not None:
            log_handle.close()
    if result.returncode != 0:
        stderr_text = result.stderr.decode("utf-8", errors="ignore") if isinstance(result.stderr, bytes) else ""
        if stderr_text:
            logging.error(stderr_text)
        if exit_on_error:
            raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd, stderr_text))
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
        logging.debug(f"Created temporary directory: {tmpdir}")
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

                if (arg_type == FullPath) or (arg_type == Optional[FullPath]):
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
            # logging.debug(f"Running {func.__name__} in temporary directory: {tmpdir}")
            os.chdir(tmpdir)

            result = func(*args, tmp_dir=tmpdir, **kwargs)
            
            # Success - clean up temp directory
            os.chdir(cwd)
            shutil.rmtree(tmpdir, ignore_errors=True)
            return result
        except Exception:
            # Error - preserve temp directory for debugging
            os.chdir(cwd)
            if uont.PRESERVE_TEMP_DIRS:
                logging.error(f"Error in {func.__name__}, preserving temp directory: {tmpdir}")
            else:
                shutil.rmtree(tmpdir, ignore_errors=True)
            raise
    return wrapper



def get_filetype(input_file):
    res = sp.run(f"htsfile {input_file}", shell=True, capture_output=True, text=True).stdout.strip()
    if "FASTQ gzip-compressed sequence data" in res:
        return "fastq.gz"
    elif "BAM version 1 compressed sequence data" in res:
        return "bam"
    else:
        raise ValueError(f"Unsupported file type for {input_file}: {res}")



g = {}



def update_dataclass(target, source):
    if type(target) is not type(source):
        raise TypeError("Both objects must be the same dataclass type")

    for field in fields(target):
        if getattr(source, field.name) is not None:
            setattr(target, field.name, getattr(source, field.name))

    return target


class JobStatus(Enum):
    NOT_RUN = "Not run"
    SUCCESS = "Success"
    FAILED = "Failed"

    def __str__(self):
        return self.value



def return_job_status(func):
    """
    Decorator to return the job status of a function.
    If the function runs successfully, it returns JobStatus.SUCCESS.
    If the function raises an exception, it returns JobStatus.FAILED.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs) -> JobStatus:
        try:
            func(*args, **kwargs)
            return JobStatus.SUCCESS
        except Exception as e:
            logging.error(f"Job {func.__name__} failed with error: {e}")
            return JobStatus.FAILED
    return wrapper

import logging
from enum import Enum
from pathlib import Path
import platform
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
import sys
from .constants import DORADO_VERSION

DEFAULT_CLI_DEPENDENCIES = [
    "autocycler",
    "bcftools",
    "chopper",
    "dnaapler",
    "dorado",
    "filtlong",
    "flye",
    "lrge",
    "medaka",
    "miniasm",
    "minimap2",
    "pigz",
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

    logging.info("Checking CLI dependencies...")
    for program in unique_programs:
        if which(program) is None:
            missing.append(program)
            sys.stdout.write(f"🔴 {program}: not found\n")
            continue

        available.append(program)
        try:
            version = get_software_version(program)
            sys.stdout.write(f"🟢 {program}: {version}\n")
        except Exception:
            sys.stdout.write(f"🟠 {program}: found, version unknown\n")

    dbs = {
        'dorado models': get_dorado_model_dir(),
        'plassembler db': get_plassembler_db_dir()
    }
    logging.info("Checking required databases...")
    for db_name, db_path in dbs.items():
        if db_path is None:
            missing.append(db_name)
            sys.stdout.write(f"🔴 {db_name}: not found\n")
        else:
            available.append(db_name)
            sys.stdout.write(f"🟢 {db_name}: {db_path}\n")

    
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


def get_dorado_model_dir() -> Optional[str]:
    dorado_executable = which("dorado")
    if dorado_executable is None:
        return None
    # check if executable is a symlink and resolve it
    dorado_executable = Path(dorado_executable)
    if dorado_executable.is_symlink():
        dorado_executable = dorado_executable.resolve()
    dorado_executable_dir = Path(os.path.dirname(dorado_executable))
    dorado_models_dir =  dorado_executable_dir.parent.parent / "dorado_models"
    if not dorado_models_dir.exists():
        return None
    return str(dorado_models_dir)


def get_plassembler_db_dir() -> Optional[str]:
    plassembler_executable = which("plassembler")
    if plassembler_executable is None:
        raise ValueError("Plassembler executable not found in PATH.")
    # check if executable is a symlink and resolve it
    plassembler_executable = Path(plassembler_executable)
    if plassembler_executable.is_symlink():
        plassembler_executable = plassembler_executable.resolve()
    plassembler_executable_dir = Path(os.path.dirname(plassembler_executable))
    plassembler_db_dir =  plassembler_executable_dir.parent / "plassembler_db"
    if not plassembler_db_dir.exists():
        return None
    return str(plassembler_db_dir)

def get_available_assemblers() -> list[str]:
    """Return a list of available assemblers from the ASSEMBLERS constant."""
    from .constants import ASSEMBLERS
    available_assemblers = []
    for assembler in ASSEMBLERS:
        if which(assembler) is not None:
            available_assemblers.append(assembler)
    return available_assemblers

def get_platform() -> str:
    """Return the platform string for Dorado download URL."""
    operating_system = {"linux":"linux","darwin":"osx"}.get(sys.platform,None)
    return operating_system

def get_architecture() -> str:
    """Return the architecture string for Dorado download URL."""
    architecture = {"x86_64":"x64","arm64":"arm64"}.get(platform.machine(),None)
    return architecture

def get_dorado_url(version: str) -> str:
    """Return the download URL for the specified version of Dorado."""
    operating_system = get_platform()
    if operating_system is None:
        raise ValueError(f"Unsupported operating system: {sys.platform}")
    architecture = get_architecture()
    if architecture is None:
        raise ValueError(f"Unsupported architecture: {platform.machine()}")


    extension = "tar.gz"  if operating_system == "linux" else "zip"
    return f"https://cdn.oxfordnanoportal.com/software/analysis/dorado-{version}-{operating_system}-{architecture}.{extension}"

def setup_dorado():
    """Check if dorado exists and if not download in the right location and add to PATH."""
    from .constants import DORADO_VERSION
    dorado_path = which("dorado")
    print(f"Checking for dorado in PATH: {dorado_path}")
    if dorado_path is None:
        # conda base dir $CONDA_PREFIX
        conda_base_dir = os.environ.get("CONDA_PREFIX")
        if conda_base_dir is None:
            raise ValueError("Dorado not found in PATH and CONDA_PREFIX is not set. Please install dorado or activate your conda environment.")
        with tempfile.TemporaryDirectory() as tmpdir:
            # download dorado from cdn
    
            url = get_dorado_url(DORADO_VERSION)
            if url.endswith(".zip"):
                cmd = f"curl -L {url} -o {tmpdir}/dorado.zip"
                sp.run(cmd, shell=True, check=True)
                # extract dorado to conda base dir
                cmd = f"unzip -q {tmpdir}/dorado.zip -d {conda_base_dir}"
                sp.run(cmd, shell=True, check=True)
            else:

                cmd = f"curl -L {url} -o {tmpdir}/dorado.tar.gz"
                sp.run(cmd, shell=True, check=True)
                # extract dorado to conda base dir
                cmd = f"tar -xzf {tmpdir}/dorado.tar.gz -C {conda_base_dir}"
                sp.run(cmd, shell=True, check=True)

            # symlink dorado to conda base dir bin
            dorado_extracted_dir = Path(conda_base_dir) / f"dorado-{DORADO_VERSION}-{get_platform()}-{get_architecture()}"
            dorado_bin_dir = dorado_extracted_dir / "bin"
            conda_bin_dir = Path(conda_base_dir) / "bin"
            symlink_path = conda_bin_dir / "dorado"
            if not symlink_path.exists():
                os.symlink(dorado_bin_dir / "dorado", symlink_path)

    if which("dorado") is None:
        raise ValueError("Dorado still not found in PATH after setup. Please check your installation.")
    else:
        if get_dorado_model_dir() is None:
            dorado_models_dir = Path(which("dorado")).parent.parent / "dorado_models"
            if not dorado_models_dir.exists():
                dorado_models_dir.mkdir(parents=True)
            cwd = os.getcwd()
            os.chdir(dorado_models_dir)
            sp.run("dorado download", shell=True, check=True)
            os.chdir(cwd)
    logging.info(f"Dorado setup complete. Dorado executable: {which('dorado')}, models directory: {get_dorado_model_dir()}")

def setup_plassembler_db():
    if get_plassembler_db_dir() is None:
        conda_base_dir = os.environ.get("CONDA_PREFIX")
        plassembler_db_dir = Path(conda_base_dir) / "plassembler_db"
        #plassembler download -d "$CONDA_PREFIX"/plassembler_db
        cmd = f"plassembler download -d {plassembler_db_dir}"
        sp.run(cmd, shell=True, check=True)

    logging.info(f"Plassembler database setup complete. Plassembler db directory: {get_plassembler_db_dir()}")
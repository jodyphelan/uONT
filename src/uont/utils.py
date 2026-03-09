import logging
import subprocess as sp
import re
from shutil import which


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
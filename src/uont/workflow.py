"""
Workflow functions for the uONT pipeline
"""

import os
import logging
from importlib import import_module
from types import SimpleNamespace
from typing import Any, Dict
from .types import FullPath

from .process import (
    process_estimate_genome_size,
    process_fastq_filter,
    process_polish,
    process_downsample_reads_to_target_depth,
    process_assemble,
    process_remove_adapters,
)


def make_dir_if_not_exists(directory: str) -> None:
    """Create a directory if it does not already exist.
    
    Args:
        directory (str): Path to the directory to create.

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def wf_assemble(
    input_fastq: FullPath,
    output_dir: str,
    tools: SimpleNamespace,
    threads: int = 4,
    min_read_depth: int = 10,
    max_contigs: int = 80,
    min_read_length: int = 1000,
    min_q_score: int = 10,
) -> None:
    """Run the assemble workflow from raw reads through polishing.

    Workflow steps:
        1. Create the output directory (if needed).
        2. Filter reads with the configured ``fastq_filter`` tool.
        3. Estimate genome size, then downsample to 150X coverage with the
           selected ``read_downsampling`` tool.
        4. Assemble the downsampled reads with the requested assembler.
        5. Polish the raw assembly using the chosen polishing tool.

    Args:
        input_fastq (FullPath): Path to the input raw fastq file.
        output_dir (str): Base output directory for workflow outputs.
        tools (SimpleNamespace): Tool selections (fastq_filter, genome_size_estimation,
            read_downsampling, assembler, polishing).
        threads (int): Number of threads to use. Defaults to 4.
        min_read_depth (int): Minimum read depth for assembly subsampling. Defaults to 10.
        max_contigs (int): Maximum number of contigs allowed in assembly. Defaults to 80.
        min_read_length (int): Minimum read length for filtering. Defaults to 1000.
        min_q_score (int): Minimum average read quality score for filtering. Defaults to 10.
    Returns:
        None
    """
    logging.info(f"Starting assembly with input {input_fastq}")
    
    make_dir_if_not_exists(f"{output_dir}/")
    
    # 1. Remove adapters
    # adapter_removed_fastq = f"{output_dir}/adapter_removed.fastq.gz"
    # process_remove_adapters(input_fastq, adapter_removed_fastq, tools.adapter_removal, threads)
    
    # 2. Filter reads by quality and length
    filtered_fastq: FullPath = FullPath(f"{output_dir}/filtered.fastq.gz")
    # process_fastq_filter(adapter_removed_fastq, filtered_fastq, threads, tools.fastq_filter)
    process_fastq_filter(input_fastq, filtered_fastq, threads, tools.fastq_filter, min_length=min_read_length, min_quality=min_q_score)
    
    # 3. Estimate genome size and Downsample to target depth
    downsampled_fastq: FullPath = FullPath(f"{output_dir}/downsampled.fastq.gz")
    genome_size = process_estimate_genome_size(filtered_fastq, tools.genome_size_estimation, threads)
    process_downsample_reads_to_target_depth(
        input_fastq=filtered_fastq, 
        output_fastq=downsampled_fastq, 
        genome_size_estimation_tool=tools.genome_size_estimation,
        read_downsampling_tool=tools.read_downsampling,
        target_depth=150, 
        genome_size=genome_size,
        threads=threads, 
    )
    
    # 4. Run assembly
    raw_assembly_file = f"{output_dir}/raw_assembly.fasta"
    process_assemble(
        input_fastq=downsampled_fastq,
        output_fasta=raw_assembly_file,
        threads=threads,
        assembler=tools.assembler,
        min_read_depth=min_read_depth,
        max_contigs=max_contigs,
        genome_size=genome_size
    )
    
    # 5. Polish assembly
    polished_assembly_file = f"{output_dir}/polished_assembly.fasta"
    process_polish(
        input_reads=filtered_fastq,
        input_assembly=raw_assembly_file,
        output_assembly=polished_assembly_file,
        threads=threads,        
        polishing_tool=tools.polishing,
    )
    
    logging.info(f"Assembly workflow completed. Final assembly: {polished_assembly_file}")


def run_configured_workflow(config: Dict[str, Any]) -> None:
    """Execute a workflow described in a YAML/JSON configuration.

    Args:
        config (dict): Parsed configuration dictionary containing a ``steps`` list.

    Raises:
        ValueError: If the configuration is malformed or references unknown steps.
    """
    steps = config.get("steps")
    if not isinstance(steps, list) or not steps:
        raise ValueError("Workflow config must define a non-empty 'steps' list")

    for idx, raw_step in enumerate(steps, start=1):
        _execute_config_step(raw_step, idx)


def _execute_config_step(raw_step: Dict[str, Any], index: int) -> None:
    if not isinstance(raw_step, dict):
        raise ValueError(f"Step {index} must be a mapping of keys like type/name/args")

    step_type = str(raw_step.get("type", "")).strip().lower()
    step_name = raw_step.get("name")
    args = raw_step.get("args") or {}

    if step_type not in {"job", "process", "workflow"}:
        raise ValueError(f"Step {index} has unsupported type '{step_type}'")
    if not step_name:
        raise ValueError(f"Step {index} must specify a 'name'")
    if not isinstance(args, dict):
        raise ValueError(f"Step {index} args must be a mapping")

    normalized_name = step_name.replace("-", "_")
    qualified_name = _qualify_step_name(step_type, normalized_name)
    func = _resolve_step_callable(step_type, qualified_name)

    logging.info(
        "[config step %s] Running %s '%s' with args=%s",
        index,
        step_type,
        qualified_name,
        args,
    )
    func(**args)


def _qualify_step_name(step_type: str, step_name: str) -> str:
    prefix_map = {"job": "job_", "process": "process_", "workflow": "wf_"}
    prefix = prefix_map[step_type]
    if step_name.startswith(prefix):
        return step_name
    return f"{prefix}{step_name}"


def _resolve_step_callable(step_type: str, qualified_name: str):
    module_name = {"job": ".jobs", "process": ".process", "workflow": ".workflow"}[step_type]
    module = import_module(module_name, package=__package__)
    if not hasattr(module, qualified_name):
        raise ValueError(f"Configured {step_type} '{qualified_name}' does not exist")
    return getattr(module, qualified_name)
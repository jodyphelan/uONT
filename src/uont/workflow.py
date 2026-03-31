"""
Workflow functions for the uONT pipeline
"""

import json
import os
import logging
from importlib import import_module
import shutil
from types import SimpleNamespace
from typing import Any, Dict

from .jobs import generate_low_dp_mask, job_collate_fasta_consensus, job_collate_flagstat_jsons, job_map_reads_minimap2, job_mapping_stats_flagstat, job_mask_low_dp_regions
from .types import FullPath

from .process import (
    process_consensus,
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
    genome_size: int = None,
    lab_id: str = None,
    link_id: str = None,
    link_directory: str = None,
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
        genome_size (int): Estimated genome size in base pairs. If not provided, will be estimated from the data.
        copy_final_assembly (str): Path to copy the final polished assembly. Defaults to None.
    Returns:
        None
    """
    logging.info(f"Starting assembly with input {input_fastq}")
    
    make_dir_if_not_exists(f"{output_dir}/")
    
    # 1. Remove adapters
    adapter_removed_fastq = f"{output_dir}/adapter_removed.fastq.gz"
    process_remove_adapters(input_fastq, adapter_removed_fastq, tools.adapter_removal, threads)
    
    # 2. Filter reads by quality and length
    filtered_fastq: FullPath = FullPath(f"{output_dir}/filtered.fastq.gz")
    process_fastq_filter(adapter_removed_fastq, filtered_fastq, threads, tools.fastq_filter, minreadlen=min_read_length, quality=min_q_score)

    
    # 3. Estimate genome size and Downsample to target depth
    downsampled_fastq: FullPath = FullPath(f"{output_dir}/downsampled.fastq.gz")
    if genome_size is None:
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

    # clean up intermediate files
    os.remove(adapter_removed_fastq)
    os.remove(filtered_fastq)
    os.remove(raw_assembly_file)
    os.remove(raw_assembly_file+".fai")
    os.remove(raw_assembly_file+".map-ont.mmi")



    if lab_id:
        final_assembly_file = f"{output_dir}/{lab_id}_polished_assembly.fasta"
        shutil.move(polished_assembly_file, final_assembly_file)
        shutil.move(downsampled_fastq, f"{output_dir}/{lab_id}_downsampled.fastq.gz")

    if link_id:
        if not os.path.exists(link_directory):
            os.makedirs(link_directory)
        os.symlink(f"{final_assembly_file}", f"{link_directory}/{link_id}_polished_assembly.fasta")

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

def wf_amplicon(
    reference_sequence: FullPath,
    input_reads: FullPath,
    output_dir: str,
    tools: SimpleNamespace,
    threads: int = 4,
    min_read_depth: int = 10,
) -> None:
    """Run a consensus generation workflow on an assembly.

    Workflow steps:
        1. Create the output directory (if needed).
        2. Polish the assembly using the selected polishing tool.

    Args:
        reference_sequence (FullPath): Path to the reference sequence fasta file.
        input_reads (FullPath): Path to the input reads fastq file.
        output_dir (str): Base output directory for workflow outputs.
        tools (SimpleNamespace): Tool selections (consensus).
        threads (int): Number of threads to use. Defaults to 4.
        min_read_depth (int): Minimum read depth for polishing subsampling. Defaults to 10.
    Returns:
        None
    """
    logging.info(f"Starting consensus generation with reference {reference_sequence} and reads {input_reads}")
    
    make_dir_if_not_exists(f"{output_dir}/")
    
    # Polish assembly
    polished_assembly_file = f"{output_dir}/polished_consensus.fasta"
    process_consensus(
        input_reads=input_reads,
        input_reference=reference_sequence,
        output_assembly=polished_assembly_file,
        threads=threads,        
        consensus_tool=tools.consensus,
    )
    bam_file = f"{output_dir}/consensus_alignment.bam"
    job_map_reads_minimap2(
        input_reads,
        polished_assembly_file,
        bam_file,
        threads=threads,
    )

    low_depth_bed_file = f"{output_dir}/low_depth_regions.bed"
    generate_low_dp_mask(
        bam_file,
        polished_assembly_file,
        low_depth_bed_file,
        min_depth=min_read_depth,
    )

    masked_consensus = f"{output_dir}/final_consensus.fasta"
    job_mask_low_dp_regions(
        input_fasta=polished_assembly_file,
        bed_file=low_depth_bed_file,
        output_fasta=masked_consensus,
    )

    flagstat_json = f"{output_dir}/flagstat.json"
    job_mapping_stats_flagstat(
        input_bam=bam_file,
        output_json=flagstat_json,
    )
    
    logging.info(f"Consensus generation workflow completed. Final consensus: {masked_consensus}")

def wf_collate_amplicon_results(
    input_directories: list[FullPath],
    output_dir: FullPath
) -> None:
    """Collate amplicon sequencing results from multiple samples into a single CSV.

    Args:
        input_directories (list[FullPath]): List of paths to input directories containing flagstat results.
        output_dir (FullPath): Path to the output directory to write collated results.

    Returns:
        None
    """
    results = []
    # make directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    job_collate_flagstat_jsons(
        input_directories=input_directories,
        output_csv=os.path.join(output_dir, "mapping_stats.csv"),
    )

    job_collate_fasta_consensus(
        input_directories=input_directories,
        output_fasta=os.path.join(output_dir, "consensus.fasta"),
    )
"""
Workflow functions for the uONT pipeline
"""

import os
import logging
from importlib import import_module
import shutil
from types import SimpleNamespace
from typing import Any, Dict, Optional

    
from .jobs import job_create_fake_asm, job_bam_to_fastq, job_dehumanise_hostile, job_nanoplot, job_ont_pre_assembly_qc, generate_low_dp_mask, job_collate_fasta_consensus, job_collate_flagstat_jsons, job_map_reads_minimap2, job_mapping_stats_flagstat, job_mask_low_dp_regions, job_remove_adapters_porechop, job_reorient_contigs_dnaapler, job_rmlst, job_write_report
from .types import FullPath

from .process import (
    process_collect_qc_results,
    process_consensus,
    process_estimate_genome_size,
    process_fastq_filter,
    process_polish,
    process_downsample_reads_to_target_depth,
    process_assemble,
    process_remove_adapters,
)

from .utils import run_in_tempdir, get_filetype, g


def make_dir_if_not_exists(directory: str) -> None:
    """Create a directory if it does not already exist.
    
    Args:
        directory (str): Path to the directory to create.

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

@run_in_tempdir
def wf_scrub(
    input_reads: FullPath,
    output_reads: FullPath,
    tools: SimpleNamespace,
    threads: int = 4,
    dehumanise: bool = True,
    sequencing_kit: Optional[str] = None,
    **kwargs
) -> None:
    """Run a QC workflow on raw reads.

    Workflow steps:
        1. Filter reads with chopper and filtlong
        2. Generate QC reports with the selected qc tool.

    Args:
        input_fastq (FullPath): Path to the input raw fastq file.
        output_fastq (FullPath): Path to the output fastq file.
        tools (SimpleNamespace): Tool selections (fastq_filter, qc).
        threads (int): Number of threads to use. Defaults to 4.
    Returns:
        None
    """
    
    logging.info(f"Starting QC workflow with input {input_reads}")
    g['workflow'] = 'qc'

    

    # 1. Filter reads
    input_filetype = get_filetype(input_reads)

    if input_filetype == "bam":
        filtered_reads = f"filtered.bam"
    else:
        filtered_reads = f"filtered.fastq.gz"

    process_remove_adapters(
        input_reads=input_reads,
        output_reads=filtered_reads,
        tool=tools.adapter_removal,
        sequencing_kit=sequencing_kit,
        threads=threads,
    )

    

    # # 2. Dehumanise reads
    # dehumanised_fastq = f"dehumanised.fastq.gz"
    # if dehumanise:
    #     job_dehumanise_hostile(
    #         input_fastq=filtered_reads,
    #         output_fastq=dehumanised_fastq,
    #         threads=threads,
    #     )
    #     final_fastq = dehumanised_fastq
    # else:
    #     final_fastq = filtered_reads

    # 3. Move selected outputs to final output directory
    shutil.move(filtered_reads, output_reads)

    logging.info(f"QC workflow completed. Filtered reads written to {output_reads}")

@run_in_tempdir
def wf_assemble(
    input_reads: FullPath,
    output_dir: FullPath,
    tools: SimpleNamespace,
    threads: int = 4,
    min_read_depth: int = 25,
    max_contigs: int = 80,
    min_read_length: int = 1000,
    min_q_score: int = 12,
    genome_size: int = None,
    rmlst: bool = False,
    lab_id: str = None,
    link_id: str = None,
    link_directory: FullPath = None,
    save_filtered_reads: bool = False,
    **kwargs
) -> None:
    """Run the assemble workflow from raw reads through polishing.

    Workflow steps:
        1. Filter reads with chopper and filtlong
        2. Estimate genome size if not provided, using the selected ``genome_size_estimation`` tool.
        3. Assemble the filtered reads with the requested assembler.
        4. Polish the raw assembly using the chosen polishing tool.
        5. Reorient the polished assembly to start at dnaA using dnaapler.

    Args:
        input_reads (FullPath): Path to the input raw reads file (bam or fastq).
        output_dir (FullPath): Base output directory for workflow outputs.
        tools (SimpleNamespace): Tool selections (fastq_filter, genome_size_estimation,
            read_downsampling, assembler, polishing).
        threads (int): Number of threads to use. Defaults to 4.
        min_read_depth (int): Minimum read depth for assembly subsampling. Defaults to 10.
        max_contigs (int): Maximum number of contigs allowed in assembly. Defaults to 80.
        min_read_length (int): Minimum read length for filtering. Defaults to 1000.
        min_q_score (int): Minimum average read quality score for filtering. Defaults to 10.
        genome_size (int): Estimated genome size in base pairs. If not provided, will be estimated from the data.
        copy_final_assembly (FullPath): Path to copy the final polished assembly. Defaults to None.
    Returns:
        None
    """
    logging.info(f"Starting assembly with input {input_reads} and output directory {output_dir}")
    make_dir_if_not_exists(f"{output_dir}/")
    
    if tools.polishing == "dorado" and not kwargs.get("bam_for_dorado"):
        raise ValueError("Dorado polishing selected but no BAM file provided. Please provide a BAM file with --bam-for-dorado.")
    
    if get_filetype(input_reads) == "bam":
        input_fastq = f"input_reads.fastq.gz"
        job_bam_to_fastq(
            input_bam=input_reads,
            output_fastq=input_fastq,
            threads=threads,
        )
    else:
        input_fastq = input_reads

        
    # 1. Filter reads
    filtered_fastq = f"filtered.fastq.gz"
    job_ont_pre_assembly_qc(
        input_fastq=input_fastq,
        output_fastq=filtered_fastq,
        threads=threads,
        quality=min_q_score,
        minreadlen=min_read_length,
    )
    
    # 2. Estimate genome size if not provided
    if genome_size is None:
        genome_size = process_estimate_genome_size(filtered_fastq, tools.genome_size_estimation, threads)
        if genome_size > 15_000_000:
            logging.error(f"Estimated genome size {genome_size} is larger than expected for a bacterial genome. Please check your data and consider providing an estimated genome size to the workflow.")
            quit()
        elif genome_size < 1_000_000:
            logging.error(f"Estimated genome size {genome_size} is smaller than expected for a bacterial genome. Please check your data and consider providing an estimated genome size to the workflow.")
            quit()
        else:
            logging.info(f"Estimated genome size: {genome_size} bp")
    
    # 3. Run assembly
    raw_assembly_file = f"raw_assembly.fasta"
    process_assemble(
        input_fastq=filtered_fastq,
        output_fasta=raw_assembly_file,
        threads=threads,
        assembler=tools.assembler,
        min_read_depth=min_read_depth,
        max_contigs=max_contigs,
        genome_size=genome_size
    )
    
    # 4. Polish assembly
    polished_assembly_file = f"polished_assembly.fasta"
    process_polish(
        input_reads=filtered_fastq,
        input_assembly=raw_assembly_file,
        output_assembly=polished_assembly_file,
        threads=threads,        
        polishing_tool=tools.polishing,
        bam_for_dorado=kwargs.get("bam_for_dorado", None),
        models_path=kwargs.get("models_directory", None),
        medaka_batch_size=kwargs.get("batch_size", None),
    )

    reoriented_assembly_file = f"polished_assembly_reoriented.fasta"
    # 5. Reorient assembly
    job_reorient_contigs_dnaapler(
        input_fasta=polished_assembly_file,
        output_fasta=reoriented_assembly_file,
        threads=threads,
    )

    job_nanoplot(
        input_fastq=input_reads,
        output_dir="nanoplot_qc",
        threads=threads,
    )
    
    # 6. write report
    run_report_file = f"run_report.json"
    job_write_report(
        input_reads=filtered_fastq,
        input_fasta=reoriented_assembly_file,
        output_report=run_report_file,
    )

    selected_outputs = {
        reoriented_assembly_file: f"{output_dir}/contigs.fasta",
        run_report_file: f"{output_dir}/run_report.json",
        "nanoplot_qc/NanoStats.txt": f"{output_dir}/NanoStats.txt",
    }
    
    if save_filtered_reads:
        selected_outputs[filtered_fastq] = f"{output_dir}/filtered_reads.fastq.gz"

    if lab_id:
        selected_outputs[reoriented_assembly_file] = f"{output_dir}/{lab_id}.contigs.fasta"
    
    for src, dst in selected_outputs.items():
        logging.info(f"Copying {src} to {dst}")
        shutil.copy(src, dst)

    if link_id:
        if not os.path.exists(link_directory):
            os.makedirs(link_directory)
        logging.info(f"Creating symlink for final assembly at {link_directory}/{link_id}.contigs.fasta")
        os.symlink(f"{reoriented_assembly_file}", f"{link_directory}/{link_id}.contigs.fasta")

    logging.info(f"Assembly workflow completed. Final assembly: {selected_outputs[reoriented_assembly_file]}")

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

@run_in_tempdir
def wf_test(
    output_dir: FullPath,
    **kwargs
):
    job_create_fake_asm(output_dir=output_dir, lab_id="test_sample")

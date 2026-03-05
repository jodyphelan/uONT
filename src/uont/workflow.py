"""
Workflow functions for the uONT pipeline
"""

import os
import logging
from types import SimpleNamespace
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
    process_fastq_filter(input_fastq, filtered_fastq, threads, tools.fastq_filter)
    
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
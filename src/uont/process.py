"""
Module containing functions that implement the core processing steps of the uONT workflow.
These functions perform specific tasks such as filtering reads, removing adapters, estimating genome size, downsampling reads, assembling genomes, and polishing assemblies. 
They are designed to be called by the higher-level workflow functions in the workflow module, which orchestrate the overall processing logic and tool selection. 
Each function takes in the necessary input files and parameters, executes the appropriate commands using the selected tools, and produces the expected output files. 
This modular design allows for flexibility in tool choice and makes it easier to maintain and extend the processing steps as needed.
"""

import os
import logging
import warnings
import sys
import pandas as pd
from dataclasses import dataclass
from tqdm import tqdm
from typing import Literal
from .utils import run_cmd
from .jobs import (
    job_assemble_autocycler,
    job_assemble_flye,
    job_consensus_bcftools,
    job_consensus_medaka,
    job_estimate_genome_size_autocycler, 
    job_fastq_filter_chopper, 
    job_remove_adapters_porechop, 
    job_downsample_filtlong,
    job_polish_medaka,
    job_estimate_genome_size_lrge,
    job_variant_calling_bcftools,
)
from .types import FullPath


@dataclass
class Sample:
    lab_id: str
    barcode: str
    run_id: str
    fastq_file: FullPath


def process_collate_barcode_fastqs(
    source_dir: str,
    output_dir: str,
    sample_sheet: str,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Collate fastq files from barcoded samples and rename according to a sample sheet.
    
    Reads an Excel sample sheet containing lab_id and barcode columns, then concatenates
    all fastq files within each barcode directory into a single file named by lab_id.
    
    Args:
        source_dir (str): Directory containing subdirectories for each barcode with fastq files.
        output_dir (str): Directory where the collated fastq files will be written.
        sample_sheet (str): Path to Excel file with lab_id and barcode columns.
    
    Returns:
        list[Sample]: Sample objects containing lab_id, barcode, and fastq_file paths.
        
    Raises:
        FileNotFoundError: If a barcode directory specified in sample sheet doesn't exist.
    """
    warnings.simplefilter(action='ignore', category=UserWarning)
    sample_sheet_df = pd.read_excel(sample_sheet, sheet_name="Sample_sheet")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        if len(os.listdir(output_dir)) > 0:
            if not force:
                logging.error(f"Output directory {output_dir} already exists and is not empty. Use --force to overwrite existing files.")
                quit(1)
            
    
    result = []
    for _, row in tqdm(list(sample_sheet_df.iterrows()), desc="Collating fastq files"):
        sample_id = row['Sample_ID']
        barcode = row['Index_barcode_#']
        run_id = row['Run_ID']
        if not os.path.exists(f'{source_dir}/{barcode}'):
            raise FileNotFoundError(f"Directory for barcode {barcode} not found in {source_dir}. Exiting.")
        fastq_files = os.listdir(f'{source_dir}/{barcode}')
        if len(fastq_files) == 0:
            logging.warning(f"No fastq files found for barcode {barcode} in {source_dir}/{barcode}. Skipping.")
            continue
        
        
        cmd = f"cat {' '.join([f'{source_dir}/{barcode}/{f}' for f in fastq_files])} > {output_dir}/{sample_id}.fastq.gz"
        if dry_run:
            for d in fastq_files:
                sys.stdout.write(f"{source_dir}/{barcode}/{d} >> {output_dir}/{sample_id}.fastq.gz\n")
        else:
            run_cmd(cmd)
        result.append(
            Sample(
                lab_id=sample_id, 
                barcode=barcode, 
                run_id=run_id,
                fastq_file=FullPath(f"{output_dir}/{sample_id}.fastq.gz")
            )
        )
    if not dry_run:
        with open(f"{output_dir}/isolates.tsv", "w") as F:
            
            for sample in result:
                F.write(f"{sample.run_id}\t{sample.lab_id}\n")
        with open(f"{output_dir}/barcodes.tsv", "w") as F:
            for sample in result:
                F.write(f"{sample.barcode}\t{sample.lab_id}\n")
    return result


def process_fastq_filter(
    input_fastq: FullPath,
    output_fastq: FullPath,
    threads: int = 4,
    tool: Literal["chopper"] = "chopper",
    quality: int = 10,
    minreadlen: int = 1000,
) -> None:
    """Filter fastq reads based on quality and length.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output filtered fastq file.
        threads (int): Number of threads to use. Defaults to 4.
        tool (Literal["chopper"]): Filtering tool to use. Defaults to "chopper".
        quality (int): Minimum quality threshold. Defaults to 10.
        minreadlen (int): Minimum read length in base pairs. Defaults to 1000.
    
    Raises:
        ValueError: If specified tool is not supported.
    """
    if tool == "chopper":
        job_fastq_filter_chopper(input_fastq, output_fastq, threads, quality, minreadlen)
    else:
        raise ValueError(f"Tool {tool} not supported for fastq filtering.")


def process_remove_adapters(
    input_fastq: FullPath,
    output_fastq: FullPath,
    tool: Literal["porechop"] = "porechop",
    threads: int = 4,
) -> None:
    """Remove sequencing adapters from reads.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output adapter-trimmed fastq file.
        tool (Literal["porechop"]): Adapter removal tool to use. Defaults to "porechop".
        threads (int): Number of threads to use. Defaults to 4.
    
    Raises:
        ValueError: If specified tool is not supported.
    """
    
    if tool == "porechop":
        job_remove_adapters_porechop(input_fastq, output_fastq, threads)
    else:
        raise ValueError(f"Tool {tool} not supported for adapter removal.")


def process_estimate_genome_size(
    input_fastq: FullPath,
    tool: Literal["autocycler"] = "autocycler",
    threads: int = 4,
) -> int:
    """Estimate genome size from sequencing reads.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        tool (Literal["autocycler"]): Tool to use for estimation. Defaults to "autocycler".
        threads (int): Number of threads to use. Defaults to 4.
    
    Returns:
        int: Estimated genome size in base pairs.
    
    Raises:
        ValueError: If specified tool is not supported.
    """
    if tool == "autocycler":
        return job_estimate_genome_size_autocycler(input_fastq, threads)
    elif tool == "lrge":
        return job_estimate_genome_size_lrge(input_fastq, threads)
    else:
        raise ValueError(f"Tool {tool} not supported for genome size estimation.")
    
def process_downsample_reads_to_target_depth(
    input_fastq: FullPath,
    output_fastq: FullPath,
    genome_size_estimation_tool: Literal["autocycler"] = "autocycler",
    read_downsampling_tool: Literal["filtlong"] = "filtlong",
    target_depth: int = 150,
    genome_size: int = None,
    threads: int = 4,
) -> None:
    """Downsample reads to achieve a target sequencing depth.

    The genome size is either supplied directly or inferred via the chosen
    ``genome_size_estimation_tool``. The target number of bases is then computed as
    ``genome_size * target_depth`` and passed to the requested
    ``read_downsampling_tool`` (currently only ``filtlong`` is supported).

    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output downsampled fastq file.
        genome_size_estimation_tool (Literal["autocycler"]): Tool used when ``genome_size`` is not supplied.
        read_downsampling_tool (Literal["filtlong"]): Implementation used to subsample reads.
        target_depth (int): Desired coverage (X). Defaults to 150.
        genome_size (int | None): Estimated genome size in base pairs, if already known.
        threads (int): Number of threads to use for estimation steps.

    Raises:
        ValueError: If the requested downsampling tool is unsupported.
    """
    if genome_size is None:
        genome_size = process_estimate_genome_size(
            input_fastq, 
            tool=genome_size_estimation_tool, 
            threads=threads
        )
    target_bases = genome_size * target_depth
    if read_downsampling_tool == "filtlong":
        job_downsample_filtlong(input_fastq, output_fastq, target_bases)
    else:
        raise ValueError(f"Tool {read_downsampling_tool} not supported for read downsampling.")
    

def process_assemble(
    input_fastq: FullPath,
    output_fasta: str,
    threads: int = 4,
    assembler: Literal["autocycler", "flye"] = "flye",
    **kwargs
) -> None:
    """Assemble reads into contigs.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fasta (str): Path to output assembly fasta file.
        threads (int): Number of threads to use. Defaults to 4.
        assembler (Literal["autocycler", "flye"]): Assembler to use. Defaults to "flye".
        **kwargs: Additional keyword arguments passed to the assembler.
    
    Raises:
        ValueError: If specified assembler is not supported.
    """
    if assembler == "autocycler":
        job_assemble_autocycler(input_fastq, output_fasta, kwargs.get("genome_size"), threads=threads, **kwargs)
    elif assembler == "flye":
        job_assemble_flye(input_fastq, output_fasta, threads, **kwargs)
    else:
        raise ValueError(f"Tool {assembler} not supported for assembly.")
    

def process_polish(
        input_reads: str,
        input_assembly: str,
        output_assembly: str,
        polishing_tool: Literal["medaka"] = "medaka",
        threads: int = 4,
) -> None:
    """Polish an assembly to improve accuracy.
    
    Uses the specified polishing tool to improve assembly quality by correcting
    errors using the original reads.
    
    Args:
        input_reads (str): Path to input fastq reads file.
        input_assembly (str): Path to input assembly fasta file to polish.
        output_assembly (str): Path where polished assembly will be written.
        polishing_tool (Literal["medaka"]): Tool to use for polishing. Defaults to "medaka".
        threads (int): Number of threads to use. Defaults to 4.
    
    Raises:
        ValueError: If specified polishing tool is not supported.
    """
    if polishing_tool == "medaka":
        job_polish_medaka(input_reads, input_assembly, output_assembly, threads)
    else:
        raise ValueError(f"Tool {polishing_tool} not supported for polishing.")


def process_consensus(
        input_reads: str,
        input_reference: str,
        output_assembly: str,
        consensus_tool: Literal["medaka","bcftools"] = "medaka",
        threads: int = 4,
) -> None:
    """Generate a consensus sequence from an assembly and reads.
    
    This function is currently identical to process_polish, but is conceptually
    distinct in that it represents the final step of the workflow where a polished
    consensus sequence is produced. In future, this function could be extended to
    include additional steps such as circularization or manual curation before polishing.
    
    Args:
        input_reads (str): Path to input fastq reads file.
        input_reference (str): Path to input assembly fasta file to generate consensus from.
        output_assembly (str): Path where consensus assembly will be written.
        consensus_tool (Literal["medaka","bcftools"]): Tool to use for polishing. Defaults to "medaka".
        threads (int): Number of threads to use. Defaults to 4.
    
    Raises:
        ValueError: If specified consensus tool is not supported.
    """
    if consensus_tool == "medaka":
        job_consensus_medaka(input_reads, input_reference, output_assembly, threads)
    elif consensus_tool == "bcftools":
        job_consensus_bcftools(input_reference, input_reads, output_assembly)
    else:
        raise ValueError(f"Tool {consensus_tool} not supported for consensus generation.")

def process_variant_calling(
    reference_fasta: str,
    input_bam: str,
    output_vcf: str,
    variant_caller: Literal["bcftools"] = "bcftools",
) -> None:
    """Call variants from mapped reads.
    
    Args:
        reference_fasta (str): Path to reference FASTA file.
        input_bam (str): Path to input BAM file with reads mapped to reference.
        output_vcf (str): Path where output VCF file will be written.
        variant_caller (Literal["bcftools"]): Tool to use for variant calling. Defaults to "bcftools".
    
    Raises:
        ValueError: If specified variant caller is not supported.
    """
    if variant_caller == "bcftools":
        job_variant_calling_bcftools(reference_fasta, input_bam, output_vcf)
    else:
        raise ValueError(f"Tool {variant_caller} not supported for variant calling.")


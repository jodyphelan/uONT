"""
All job functions that perform specific tasks in the uONT workflow are defined here. 
Each function corresponds to a specific step in the processing pipeline, 
such as filtering reads, removing adapters, assembling genomes, 
polishing assemblies, etc. These functions are designed to be called by the 
higher-level process functions in the process module, which handle the overall 
workflow logic and tool selection. They are also served as command implementations 
for the CLI interface, allowing users to run specific steps directly from the 
command line if desired.
"""

import os
import logging
import tempfile
import shutil
from pathogenprofiler_core import run_cmd
import functools
from .types import FullPath

def run_in_tempdir(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        try:
            os.chdir(tmpdir)
            return func(*args, tmp_dir=tmpdir, **kwargs)
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir, ignore_errors=True)
    return wrapper

def job_fastq_filter_chopper(
    input_fastq: FullPath,
    output_fastq: FullPath,
    threads: int = 4,
    quality: int = 10,
    minreadlen: int = 1000,
) -> None:
    """Filter and quality-trim reads using chopper.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output filtered fastq file.
        threads (int): Number of threads to use. Defaults to 4.
        quality (int): Minimum quality threshold. Defaults to 10.
        minreadlen (int): Minimum read length in base pairs. Defaults to 1000.

    Returns:
        None
    """
    logging.info(f"Running chopper on {input_fastq} with output {output_fastq} using {threads} threads.")
    cmd = f"chopper -t {threads} -q {quality} -l {minreadlen} -i {input_fastq} > {output_fastq}"
    run_cmd(cmd)

# src/uont/jobs.py
@run_in_tempdir
def job_assemble_miniasm(
    input_fastq: FullPath, 
    output_assembly: FullPath, 
    threads: int,
    **kwargs
) -> None:
    """
    Assemble reads using miniasm.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_assembly: Path where the assembly FASTA will be written.
        genome_size: Estimated genome size for assembly parameters.
        threads: Number of threads to use for assembly.
    """
    logging.info(f"Running miniasm with input: {input_fastq}")

    # First we run minimap2 to generate the overlaps (not shown here), then we run miniasm on the overlaps.
    cmd = f"minimap2 -x ava-ont -t {threads} {input_fastq} {input_fastq} | pigz -c > overlaps.paf.gz"
    run_cmd(cmd)

    # Construct and log the miniasm command here
    cmd = f"miniasm -f {input_fastq} overlaps.paf.gz > miniasm_assembly.gfa"
    # Execute the command (e.g., via subprocess)
    run_cmd(cmd)

    # The output of miniasm is a GFA file. We need to convert it to FASTA format. 
    for line in open("miniasm_assembly.gfa"):
        if line.startswith("S"):
            parts = line.strip().split("\t")
            contig_name = parts[1]
            contig_seq = parts[2]
            with open("miniasm_assembly.fasta", "a") as fasta_out:
                fasta_out.write(f">{contig_name}\n{contig_seq}\n")
     
    # Copy the output assembly to the final destination if needed (e.g., if miniasm writes to a temp file)
    shutil.copy("miniasm_assembly.fasta", output_assembly)

def job_remove_adapters_porechop(
    input_fastq: FullPath,
    output_fastq: FullPath,
    threads: int = 4,
) -> None:
    """Remove sequencing adapters using porechop_abi.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output adapter-trimmed fastq file.
        threads (int): Number of threads to use. Defaults to 4.

    Returns:
        None
    """
    logging.info(f"Running porechop on {input_fastq} with output {output_fastq} using {threads} threads.")
    cmd = f"porechop_abi -t {threads} -i {input_fastq} -o {output_fastq}"
    run_cmd(cmd)

@run_in_tempdir
def job_assemble_flye(
    input_fastq: FullPath,
    output_fasta: FullPath,
    threads: int = 4,
    **kwargs
) -> None:
    """Assemble reads into contigs using flye.
    
    Args:
        input_fastq (FullPath): Path to input fastq file with reads.
        output_fasta (FullPath): Path to output assembly fasta file.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Optional Flye parameters such as ``genome_size``.

    Returns:
        None
    """
    
    logging.info(f"Running flye assembly on {input_fastq} using {threads} threads.")
        
    genome_size_string = f"--genome-size {kwargs['genome_size']}" if 'genome_size' in kwargs else ""
    cmd = f"flye --nano-raw {input_fastq} --out-dir assembly --threads {threads} {genome_size_string}"
    run_cmd(cmd)
    
    # move final assembly to output location
    shutil.move(f"assembly/assembly.fasta", output_fasta)

@run_in_tempdir
def job_assemble_autocycler(
    input_fastq: FullPath,
    output_fasta: FullPath,
    genome_size: int,
    threads: int = 4,
    assembler: str = "miniasm",
    min_read_depth: int = 10,
    max_contigs: int = 80,
    **kwargs
) -> None:
    """Run complete autocycler assembly workflow.
    
    Performs the full autocycler pipeline including subsampling, multiple assemblies,
    compression, clustering, trimming, resolving, and combining into final consensus.
    
    Workflow:
        1. Subsample reads into multiple subsets
        2. Run QC on subsampled reads with seqkit
        3. Assemble each subset using specified assembler
        4. Compress assemblies
        5. Cluster assemblies
        6. Trim and resolve clusters
        7. Combine clusters into single consensus assembly
    
    Args:
        input_fastq (FullPath): Path to input fastq file with reads.
        output_fasta (FullPath): Destination path for the consensus FASTA file.
        genome_size (int): Estimated genome size in base pairs.
        threads (int): Number of threads to use. Defaults to 4.
        assembler (str): Assembler helper to use (``miniasm``, ``flye``, ``raven``).
        min_read_depth (int): Minimum read depth for subsampling. Defaults to 10.
        max_contigs (int): Maximum number of contigs allowed. Defaults to 80.
        kwargs (dict[str, object]): Extra parameters accepted for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Starting autocycler assembly workflow")
    
   
    # Create output directory structure
    autocycler_output_dir = 'autocycler'
    if not os.path.exists(autocycler_output_dir):
        os.makedirs(autocycler_output_dir)
    
    # 1. Subsample reads
    logging.info(f"Subsampling reads into subsets with minimum read depth {min_read_depth}X")
    subsample_cmd = f"autocycler subsample --reads {input_fastq} --out_dir {autocycler_output_dir}/subsampled_reads --genome_size {genome_size} --min_read_depth {min_read_depth}"
    run_cmd(subsample_cmd)
    
    # 2. QC on subsampled reads
    logging.info("Running QC on subsampled reads")
    qc_cmd = f"seqkit stats -Ta {autocycler_output_dir}/subsampled_reads/sample_0*.fastq > {autocycler_output_dir}/subsamples_QC.tsv"
    run_cmd(qc_cmd)
    
    # 3. Make assemblies for each subsampled read
    logging.info(f"Creating assemblies using {assembler}")
    os.makedirs(f"{autocycler_output_dir}/autocycler_assemblies", exist_ok=True)
    
    # Get list of sample files
    import glob
    sample_files = sorted(glob.glob(f"{autocycler_output_dir}/subsampled_reads/sample_*.fastq"))
    
    for sample_file in sample_files:
        sample_num = os.path.basename(sample_file).replace("sample_", "").replace(".fastq", "")
        logging.info(f"Assembling sample {sample_num} with {assembler}")
        assembly_cmd = f"autocycler helper {assembler} --reads {sample_file} --out_prefix {autocycler_output_dir}/autocycler_assemblies/{assembler}_{sample_num} --threads {threads} --genome_size {genome_size}"
        run_cmd(assembly_cmd)
    
    # 4. Compress assemblies
    logging.info("Compressing assemblies")
    compress_cmd = f"autocycler compress -i {autocycler_output_dir}/autocycler_assemblies -a {autocycler_output_dir}/autocycler_out --threads {threads} --max_contigs {max_contigs}"
    run_cmd(compress_cmd)
    
    # 5. Cluster assemblies
    logging.info("Clustering assemblies")
    cluster_cmd = f"autocycler cluster -a {autocycler_output_dir}/autocycler_out --max_contigs {max_contigs}"
    run_cmd(cluster_cmd)
    
    # 6. Trim and resolve clusters
    logging.info("Trimming and resolving clusters")

    ## check for qc_pass clusters
    if not glob.glob(f"{autocycler_output_dir}/autocycler_out/clustering/qc_pass/cluster_*"):
        logging.warning("No clusters passed QC. Check the clustering output and adjust parameters if necessary.")
        return None
    
    cluster_dirs = glob.glob(f"{autocycler_output_dir}/autocycler_out/clustering/qc_pass/cluster_*")
    for cluster_dir in sorted(cluster_dirs):
        logging.info(f"Processing {os.path.basename(cluster_dir)}")
        trim_cmd = f"autocycler trim -c {cluster_dir} --threads {threads}"
        run_cmd(trim_cmd)
        resolve_cmd = f"autocycler resolve -c {cluster_dir}"
        run_cmd(resolve_cmd)
    
    # 7. Combine clusters into single assembly
    logging.info("Combining clusters into single assembly")
    combine_cmd = f"autocycler combine -a {autocycler_output_dir}/autocycler_out -i {autocycler_output_dir}/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa"
    run_cmd(combine_cmd)

    # Move final assembly to output location
    shutil.move(f"{autocycler_output_dir}/autocycler_out/consensus_assembly.fasta", output_fasta)
    logging.info(f"Autocycler assembly workflow completed. Final assembly written to {output_fasta}")

def job_estimate_genome_size_autocycler(
    input_fastq: FullPath,
    threads: int = 4,
) -> int:
    """Estimate genome size using autocycler.
    
    Uses autocycler's helper genome_size command to estimate the genome size
    from the input reads.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        threads (int): Number of threads to use. Defaults to 4.
    
    Returns:
        int: Estimated genome size in base pairs.
    """
    logging.info(f"Running autocycler genome size estimation on {input_fastq} using {threads} threads.")
    #autocycler helper genome_size --reads $isolate.trim.q"$q_out".1000.fastq --threads "$t"
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        cmd = f"autocycler helper genome_size --reads {input_fastq} --threads {threads} > {tmp.name}"
        run_cmd(cmd)
        genome_size = tmp.read().decode().strip()
        logging.info(f"Estimated genome size: {genome_size}")
    return int(genome_size)


def job_downsample_filtlong(
    input_fastq: FullPath,
    output_fastq: FullPath,
    target_bases: int,
) -> None:
    """Downsample reads to a target number of bases using filtlong.
    
    Selects the highest quality reads to reach the target base count,
    using mean quality weighting.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output downsampled fastq file.
        target_bases (int): Target number of bases to retain.

    Returns:
        None
    """
    logging.info(f"Running filtlong downsample on {input_fastq} with target bases {target_bases} to output {output_fastq}.")
    cmd = f"filtlong --target_bases {target_bases} --mean_q_weight 10 {input_fastq} | pigz -c > {output_fastq}"
    run_cmd(cmd)


@run_in_tempdir
def job_polish_medaka(
        input_reads: FullPath,
        input_assembly: FullPath,
        output_assembly: FullPath,
        threads: int = 4,
        **kwargs
) -> None:
    """Polish an assembly twice using medaka consensus.
    
    Performs two rounds of medaka polishing with bacterial model to improve
    assembly accuracy. The first round polishes the input assembly, and the
    second round polishes the first round output.
    
    Args:
        input_reads (FullPath): Path to input fastq reads file used for polishing.
        input_assembly (FullPath): Path to input assembly fasta file to polish.
        output_assembly (FullPath): Path where final polished assembly will be written.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Additional options including ``tmp_dir``.

    Returns:
        None
    """
    
    logging.info(f"Running medaka polishing (round 1) for {input_assembly}")
    
    tmp_dir = kwargs.get("tmp_dir")
    tmp_dir1 = os.path.join(tmp_dir, "round1")
    os.makedirs(tmp_dir1, exist_ok=True)
    # Round 1: Polish original assembly
    cmd_round1 = f"medaka_consensus -t {threads} -i {input_reads} -d {input_assembly} -o {tmp_dir1} --bacteria"
    run_cmd(cmd_round1)

    tmp_dir2 = os.path.join(tmp_dir, "round2")
    os.makedirs(tmp_dir2, exist_ok=True)
    cmd_round2 = f"medaka_consensus -t {threads} -i {input_reads} -d {tmp_dir1}/consensus.fasta -o {tmp_dir2} --bacteria"
    run_cmd(cmd_round2)

    shutil.move(f"{tmp_dir2}/consensus.fasta", output_assembly)
    logging.info(f"Medaka polishing completed. Output: {output_assembly}")

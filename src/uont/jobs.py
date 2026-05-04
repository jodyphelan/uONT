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

import json
import os
import logging
import tempfile
import shutil
import pysam
import numpy as np
from tqdm import tqdm
import platform
import subprocess as sp
from typing import Tuple
from .utils import run_cmd, run_in_tempdir, timeit
from .types import FullPath
from .qc import Fasta




# def run_in_tempdir(func):
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         cwd = os.getcwd()
#         tmpdir = tempfile.mkdtemp()
#         try:
#             # filter out any arguments in kwargs that are already in args to avoid duplication
#             # move kwargs to args based on function signature
#             sig = inspect.signature(func)
#             # find positional arguments in kwargs and move them to args

#             if len(args)==0:
#                 new_args = {}
#                 new_kwargs = {}
#                 for param in sig.parameters.values():
#                     if param.name=='kwargs':
#                         continue
#                     if param.default is inspect.Parameter.empty:
#                         new_args[param.name] = kwargs.pop(param.name)
#                     else:
#                         new_kwargs[param.name] = kwargs.get(param.name, param.default)
#                 args = tuple(new_args.values())
#                 kwargs = new_kwargs
#             else:
#                 new_args = args
#                 new_kwargs = kwargs
            
#             kwargs = {k: v for k, v in kwargs.items() if k not in func.__code__.co_varnames}
#             # convert any FullPath arguments to absolute paths
#             sig = inspect.signature(func)
#             for param in sig.parameters.values():
#                 arg_type = param.annotation if param.annotation is not param.empty else str
#                 if arg_type == FullPath:
#                     arg_index = list(sig.parameters).index(param.name)
#                     if arg_index < len(args):
#                         args = list(args)
#                         args[arg_index] = os.path.abspath(args[arg_index])
#                         args = tuple(args)
#                     elif param.name in kwargs:
#                         kwargs[param.name] = os.path.abspath(kwargs[param.name])

#             os.chdir(tmpdir)
#             return func(*args, tmp_dir=tmpdir, **kwargs)
#         finally:
#             os.chdir(cwd)
#             shutil.rmtree(tmpdir, ignore_errors=True)
#     return wrapper

@timeit
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
    cmd = f"chopper -t {threads} -q {quality} -l {minreadlen} -i {input_fastq} | pigz -p {threads} -c > {output_fastq}"
    run_cmd(cmd)

# src/uont/jobs.py

@timeit
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
    cmd = f"minimap2 -x ava-ont -t {threads} {input_fastq} {input_fastq} | pigz -p {threads} -c > overlaps.paf.gz"
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

def job_ont_pre_assembly_qc(
    input_fastq: FullPath,
    output_fastq: FullPath,
    threads: int = 4,
    quality: int = 12,
    minreadlen: int = 1000,
    headcrop: int = 80,
    tailcrop: int = 80,
    keeppercent: int = 90,
) -> None:
    """Basic filtering of nanopore reads using chopper and filtlong.
    
    Args:
        input_fastq (FullPath): Path to input fastq file.
        output_fastq (FullPath): Path to output filtered fastq file.
        threads (int): Number of threads to use. Defaults to 4.
        quality (int): Minimum quality threshold. Defaults to 10.
        minreadlen (int): Minimum read length in base pairs. Defaults to 1000.

    Returns:
        None
    """
    logging.info(f"Running basic filtering on {input_fastq} with output {output_fastq} using {threads} threads.")
    # chopper
    cmd = f"gunzip -c {input_fastq} | chopper --minlength {minreadlen} --quality {quality} --headcrop {headcrop} --tailcrop {tailcrop} -t {threads} > intermediate.fastq"
    run_cmd(cmd)
    # filtlong
    cmd = f"filtlong --min_length {minreadlen} --keep_percent {keeppercent} intermediate.fastq | pigz -p {threads} -c > {output_fastq}"
    run_cmd(cmd)

@timeit
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

@timeit
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
def job_reorient_contigs_dnaapler(
    input_fasta: FullPath,
    output_fasta: FullPath,
    threads: int = 4,
    **kwargs
):
    """Reorient contigs to start at dnaA using dnaapler.
    
    Args:
        input_fasta (FullPath): Path to input assembly fasta file.
        output_fasta (FullPath): Path where reoriented assembly fasta will be written.
        threads (int): Number of threads to use. Defaults to 4.

    Returns:
        None
    """
    logging.info(f"Running dnaapler to reorient contigs in {input_fasta} using {threads} threads.")
    cmd = f"dnaapler all --input {input_fasta} --output output.dnaapler --threads {threads}"
    run_cmd(cmd)
    shutil.move("output.dnaapler/dnaapler_reoriented.fasta", output_fasta)

@timeit
@run_in_tempdir
def job_assemble_autocycler(
    input_fastq: FullPath,
    output_fasta: FullPath,
    genome_size: int,
    threads: int = 4,
    assemblers: Tuple[str] = ("flye", "miniasm","nextdenovo", "raven"),
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
        assemblers (Tuple[str]): Assemblers to use (``miniasm``, ``flye``, ``raven``). Defaults to ``("flye", "miniasm")``.
        min_read_depth (int): Minimum read depth for subsampling. Defaults to 10.
        max_contigs (int): Maximum number of contigs allowed. Defaults to 80.
        kwargs (dict[str, object]): Extra parameters accepted for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Starting autocycler assembly workflow")
    
    if platform.machine() == "arm64":
        assemblers = tuple(a for a in assemblers if a != "nextdenovo")
        logging.warning("NextDenovo is not currently supported on ARM64 architecture. It will be skipped in the autocycler workflow.")
   
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
    logging.info(f"Creating assemblies using {assemblers}")
    os.makedirs(f"{autocycler_output_dir}/autocycler_assemblies", exist_ok=True)
    
    # Get list of sample files
    import glob
    sample_files = sorted(glob.glob(f"{autocycler_output_dir}/subsampled_reads/sample_*.fastq"))
    
    combinations = [(sample_file, assembler) for sample_file in sample_files for assembler in assemblers]
    for sample_file, assembler in tqdm(combinations, desc="Assembling subsamples"):
        sample_num = os.path.basename(sample_file).replace("sample_", "").replace(".fastq", "")
        logging.debug(f"Assembling sample {sample_num} with {assembler}")
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
    for cluster_dir in tqdm(sorted(cluster_dirs),desc="Processing clusters"):
        logging.debug(f"Processing {os.path.basename(cluster_dir)}")
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

@timeit
def job_estimate_genome_size_lrge(
    input_fastq: FullPath,
    threads: int = 4,
):
    """Estimate genome size using lrge.

    Uses lrge's estimate_genome_size command to estimate the genome size
    from the input reads.

    Args:
        input_fastq (FullPath): Path to input fastq file.
        threads (int): Number of threads to use. Defaults to 4.

    Returns:
        int: Estimated genome size in base pairs.
    """
    logging.info(f"Running lrge genome size estimation on {input_fastq} using {threads} threads.")
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        cmd = f"lrge {input_fastq} -t {threads} -o {tmp.name}"
        run_cmd(cmd)
        genome_size = tmp.read().decode().strip()
        logging.info(f"Estimated genome size: {genome_size}")
    return int(genome_size)

@timeit
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


@timeit
def job_downsample_filtlong(
    input_fastq: FullPath,
    output_fastq: FullPath,
    target_bases: int,
    threads: int = 4,
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
    cmd = f"filtlong --target_bases {target_bases} --mean_q_weight 10 {input_fastq} | pigz -p {threads} -c > {output_fastq}"
    run_cmd(cmd)

@timeit
@run_in_tempdir
def job_polish_dorado(
    input_bam: FullPath,
    input_assembly: FullPath,
    output_assembly: FullPath,
    threads: int = 4,
    **kwargs
) -> None:
    """Polish an assembly using dorado.
    
    Args:
        input_bam (FullPath): Path to input BAM file used for polishing.
        input_assembly (FullPath): Path to input assembly fasta file to polish.
        output_assembly (FullPath): Path where polished assembly will be written.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Additional options including ``tmp_dir``.

    Returns:
        None
    """
    logging.info(f"Running dorado polishing on {input_assembly}")
    
    input_bam = kwargs.get("bam_for_dorado", input_bam)

    tmp_dir = kwargs.get("tmp_dir")
    os.makedirs(tmp_dir, exist_ok=True)
    
    cmd = f"dorado  aligner -t 4 {input_assembly} {input_bam} | samtools sort -@ 4 -o aligned.bam -"
    run_cmd(cmd)

    cmd = f"samtools index aligned.bam"
    run_cmd(cmd)

    cmd = f"dorado  polish -t {threads} --ignore-read-groups --bacteria aligned.bam {input_assembly} > {output_assembly}"
    run_cmd(cmd)



@timeit
@run_in_tempdir
def job_polish_medaka(
    input_reads: FullPath,
    input_assembly: FullPath,
    output_assembly: FullPath,
    threads: int = 4,
    batch_size: int = 100,
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
    cmd_round1 = f"OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 medaka_consensus -b {batch_size} -t {threads} -i {input_reads} -d {input_assembly} -o {tmp_dir1} --bacteria"
    # run command and check check stdout for "ERROR"
    result1 = sp.run(cmd_round1, shell=True, capture_output=True, text=True)
    bacteria_model_working = True
    if result1.returncode != 0:
        if "ERROR: --bacteria was specified but input data was not compatible." in result1.stdout:
            logging.warning("Medaka round 1 failed with --bacteria model. Retrying without --bacteria flag.")
            bacteria_model_working = False
            cmd_round1 = f"OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 medaka_consensus -b {batch_size} -t {threads} -i {input_reads} -d {input_assembly} -o {tmp_dir1}"
            run_cmd(cmd_round1)
        else:
            logging.error(f"Medaka round 1 failed with error:\n{result1.stdout}")
            raise ValueError("Medaka polishing failed in round 1.")

    tmp_dir2 = os.path.join(tmp_dir, "round2")
    os.makedirs(tmp_dir2, exist_ok=True)
    if bacteria_model_working:
        logging.info(f"Running medaka polishing (round 2) with bacteria model for {tmp_dir1}/consensus.fasta")
        cmd_round2 = f"OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 medaka_consensus -b {batch_size} -t {threads} -i {input_reads} -d {tmp_dir1}/consensus.fasta -o {tmp_dir2} --bacteria"
    else:
        cmd_round2 = f"OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 medaka_consensus -b {batch_size} -t {threads} -i {input_reads} -d {tmp_dir1}/consensus.fasta -o {tmp_dir2}"
    run_cmd(cmd_round2)

    shutil.move(f"{tmp_dir2}/consensus.fasta", output_assembly)
    logging.info(f"Medaka polishing completed. Output: {output_assembly}")

@run_in_tempdir
def job_map_reads_minimap2(
    input_reads: FullPath,
    input_assembly: FullPath,
    output_bam: FullPath,
    threads: int = 4,
    **kwargs
) -> None:
    """Map reads to an assembly using minimap2 and output sorted BAM.
    
    Args:
        input_reads (FullPath): Path to input fastq reads file.
        input_assembly (FullPath): Path to input assembly fasta file.
        output_bam (FullPath): Path where sorted BAM file will be written.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Mapping reads to assembly with minimap2. Reads: {input_reads}, Assembly: {input_assembly}, Output BAM: {output_bam}")
    cmd = f"minimap2 -ax map-ont -t {threads} {input_assembly} {input_reads} | samtools sort -o {output_bam}"
    run_cmd(cmd)
    cmd = f"samtools index {output_bam}"
    run_cmd(cmd)

def job_mask_low_dp_regions(
    input_fasta: FullPath,
    bed_file: FullPath,
    output_fasta: FullPath
) -> None:
    """Mask low depth regions in an assembly based on a BED file.
    
    Replaces sequences in the input FASTA with Ns at positions specified in the BED file.
    
    Args:
        input_fasta (FullPath): Path to input assembly fasta file.
        bed_file (FullPath): Path to BED file containing regions to mask.
        output_fasta (FullPath): Path where masked assembly fasta will be written.

    Returns:
        None
    """
    logging.info(f"Masking low depth regions in {input_fasta} using BED file {bed_file}. Output: {output_fasta}")
    refseq = pysam.FastaFile(input_fasta)
    seqname = refseq.references[0]
    sequence = list(refseq.fetch(seqname))
    
    for line in open(bed_file):
        row = line.strip().split("\t")
        start, end = int(row[1]), int(row[2])
        for i in range(start, end):
            sequence[i] = "N"
    
    with open(output_fasta, "w") as O:
        O.write(f">{seqname}\n{''.join(sequence)}\n")


@run_in_tempdir
def job_rmlst(
    input_fasta: FullPath,
    output_tsv: FullPath,
    **kwargs
) -> None:
    """Run rMLST analysis on an assembly to determine sequence type.
    
    Args:
        input_fasta (FullPath): Path to input assembly fasta file.
        output_tsv (FullPath): Path where rMLST results will be written.

    Returns:
        None
    """
    logging.info(f"Running rMLST analysis on {input_fasta}. Output: {output_tsv}")
    cmd = f"rmlst -f {input_fasta} --species-only -o {output_tsv}"
    run_cmd(cmd)

@run_in_tempdir
def generate_low_dp_mask(
    bam: FullPath,
    ref: FullPath,
    outfile: FullPath,
    min_dp: int = 10,
    **kwargs
) -> None:
    """
    Generate a BED file of low depth regions based on a BAM file and reference FASTA.
    Args:
        bam (FullPath): Path to input BAM file with reads mapped to reference.
        ref (FullPath): Path to reference FASTA file.
        outfile (FullPath): Path where output BED file will be written.
        min_dp (int): Minimum read depth threshold for masking. Defaults to 10.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    refseq = pysam.FastaFile(ref)
    seqname = refseq.references[0]
    dp = np.zeros(refseq.get_reference_length(seqname))
    tmp_depth_file = f"{outfile}.depth"
    cmd = f"samtools depth {bam} > {tmp_depth_file}"
    run_cmd(cmd)
    for l in open(tmp_depth_file):
        row = l.strip().split("\t")
        dp[int(row[1])-1] = int(row[2])

    masked_positions = []
    for i, d in enumerate(dp):
        if d < min_dp:
            masked_positions.append(i)

    with open(outfile,"w") as O:
        for i in masked_positions:
            O.write(f"{seqname}\t{i}\t{i+1}\n")

    shutil.move(outfile, f"{outfile}")

@timeit
@run_in_tempdir
def job_dehumanise_hostile(
        input_fastq: FullPath,
        output_fastq: FullPath,
        threads: int = 4,
        **kwargs
) -> None:
    """Dehumanise reads using hostile.
    
    Args:
        input_fastq (FullPath): Path to input fastq file with reads.
        output_fastq (FullPath): Path where dehumanised fastq will be written.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Running hostile to dehumanise reads in {input_fastq}. Output: {output_fastq}")
    gzip_output = output_fastq.endswith(".gz")
    if gzip_output:
        output_fastq = output_fastq[:-3]
    cmd = f"hostile clean --fastq1 {input_fastq} -t {threads} -o - > {output_fastq}"
    run_cmd(cmd)
    if gzip_output:
        cmd = f"pigz -p {threads} {output_fastq}"
        run_cmd(cmd)


def job_qc_python(
    input_fasta: FullPath,
    output_tsv: FullPath,
    sample_id: str = None,
    **kwargs
) -> None:
    """Run QC on reads using a Python-based approach.
    
    This is a placeholder for a custom QC implementation that could be developed
    in Python, potentially using libraries like Biopython or SeqIO to analyze
    read quality and length distributions, and output a QC report.

    Args:
        input_fasta (FullPath): Path to input fasta file with reads.
        output_tsv (FullPath): Path where QC report will be written.
        sample_id (str): Sample identifier. Defaults to None.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Running custom Python QC on {input_fasta}. Output: {output_tsv}")
    fasta = Fasta(input_fasta, sample_id=sample_id)
    fasta.write_qc_report(output_file=output_tsv)



def job_variant_calling_bcftools(
    reference_fasta: FullPath,
    input_bam: FullPath,
    output_vcf: FullPath,
    **kwargs
) -> None:
    """Generate a VCF file using bcftools.
    
    Args:
        reference_fasta (FullPath): Path to reference FASTA file.
        input_bam (FullPath): Path to input BAM file with reads mapped to reference.
        output_vcf (FullPath): Path where output VCF file will be written.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    logging.info(f"Generating VCF with bcftools. BAM: {input_bam}, Reference: {reference_fasta}, Output VCF: {output_vcf}")
    cmd = f"bcftools mpileup -a AD -f {reference_fasta} {input_bam} | bcftools call -mv -Oz -o {output_vcf}"
    run_cmd(cmd)
    cmd = f"bcftools index {output_vcf}"
    run_cmd(cmd)


def job_consensus_bcftools(
    reference_fasta: FullPath,
    input_fastq: FullPath,
    output_fasta: FullPath,
    **kwargs
) -> None:
    """Generate a consensus FASTA from a VCF file using bcftools.
    
    Args:
        reference_fasta (FullPath): Path to reference FASTA file.
        input_vcf (FullPath): Path to input VCF file with variants.
        output_fasta (FullPath): Path where output consensus FASTA will be written.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    bam_filename = f"{output_fasta}.bam"
    job_map_reads_minimap2(
        input_reads=input_fastq,
        input_assembly=reference_fasta,
        output_bam=bam_filename,
        threads=kwargs.get("threads", 4)
    )
    vcf_filename = f"{output_fasta}.vcf.gz"
    job_variant_calling_bcftools(
        reference_fasta=reference_fasta,
        input_bam=bam_filename,
        output_vcf=vcf_filename,
    )
    annotated_vcf_filename = f"{output_fasta}.annotated.vcf.gz"
    job_vcf_annotate_maaf(
        input_vcf=vcf_filename,
        output_vcf=annotated_vcf_filename
    )
    filtered_vcf_filename = f"{output_fasta}.filtered.vcf.gz"
    cmd = f"bcftools view -i 'INFO/MAAF>0.5' {annotated_vcf_filename} -Oz -o {filtered_vcf_filename}"
    run_cmd(cmd)
    cmd = f"bcftools index {filtered_vcf_filename}"
    run_cmd(cmd)
    logging.info(f"Generating consensus FASTA with bcftools. VCF: {filtered_vcf_filename}, Reference: {reference_fasta}, Output FASTA: {output_fasta}")
    cmd = f"bcftools consensus -f {reference_fasta} {filtered_vcf_filename} > {output_fasta}"
    run_cmd(cmd)

def job_consensus_medaka(
    input_reads: FullPath,
    input_assembly: FullPath,
    output_assembly: FullPath,
    threads: int = 4,
    **kwargs
):
    """Generate a consensus FASTA from an assembly and reads using medaka.
    
    Args:
        input_reads (FullPath): Path to input fastq reads file.
        input_assembly (FullPath): Path to input assembly fasta file to generate consensus from.
        output_assembly (FullPath): Path where consensus assembly will be written.
        threads (int): Number of threads to use. Defaults to 4.
        kwargs (dict[str, object]): Additional options for interface symmetry.

    Returns:
        None
    """
    job_polish_medaka(
        input_reads=input_reads,
        input_assembly=input_assembly,
        output_assembly=output_assembly,
        threads=threads,
        **kwargs
    )

def job_vcf_annotate_maaf(
    input_vcf: str,
    output_vcf: str
):

    vcf_in = pysam.VariantFile(input_vcf, "r")
    header = vcf_in.header

    # add float INFO field
    tag = {'ID':'MAAF','Number':'1','Type':'Float','Description':'Major Alternative Allele Frequency'}
    header.add_meta('INFO',items=tag.items())


    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)
    for var in vcf_in:
        # short variant
        if 'AD' in var.samples[0]:
            if var.samples[0]['AD'] == (None,):
                maaf = 0
            else:
            # calculate MAF from FMT/AD
                maaf = max(var.samples[0]['AD'][1:])/sum(var.samples[0]['AD'])

        elif 'SVTYPE' in var.info:
            # long variant
            # use DR:DV:RR:RV
            total_reads = var.samples[0]['DR'] + var.samples[0]['DV'] + var.samples[0]['RR'] + var.samples[0]['RV']
            if var.samples[0]['RR'] == (None,):
                maaf = 0
            elif total_reads==0:
                maaf = 0
            else:
                maaf = (var.samples[0]['DV'] + var.samples[0]['RV'])/total_reads
        else:
            # panic
            raise ValueError("No AD or SVTYPE in INFO field")
        var.info['MAAF'] = maaf
        vcf_out.write(var)

    vcf_in.close()
    vcf_out.close()

def job_mapping_stats_flagstat(
    input_bam: FullPath,
    output_json: FullPath
):
    """Generate mapping statistics in JSON format using samtools flagstat.
    
    Args:
        input_bam (FullPath): Path to input BAM file with reads mapped to reference.
        output_json (FullPath): Path where output JSON file with mapping stats will be written.
    Returns:
        None
    """
    logging.info(f"Generating mapping statistics with samtools flagstat. Input BAM: {input_bam}, Output JSON: {output_json}")
    cmd = f"samtools flagstat {input_bam} -O json > {output_json}"
    run_cmd(cmd)

def extract_flagstat_stats(
    input_json: FullPath
) -> dict:
    """Extract relevant statistics from a samtools flagstat JSON output.
    
    Args:
        input_json (FullPath): Path to input JSON file generated by samtools flagstat.
    Returns:
        dict: Dictionary containing extracted statistics such as total reads, mapped reads, and percentage mapped.
    """
    with open(input_json) as f:
        data = json.load(f)
        total = data["QC-passed reads"]["total"]
        mapped = data["QC-passed reads"]["mapped"]
        percent_mapped = mapped/total*100 if total > 0 else 0
        return {
            'total': total,
            'mapped': mapped,
            '% mapped': percent_mapped
        }

def job_collate_flagstat_jsons(
    input_directories: list[FullPath],
    output_csv: FullPath
):

# extract the "QC-passed reads": ["total","mapped"]  from each json and average them across all jsons, then write to output_json
    results = []
    for dir in input_directories:
        json_path = os.path.join(dir, "flagstat.json")
        stats = extract_flagstat_stats(json_path)
        sample_name = os.path.basename(dir)
        results.append({
            'sample_id': sample_name,
            'total': stats['total'],
            'mapped': stats['mapped'],
            '% mapped': stats['% mapped']
        })

    with open(output_csv, "w") as O:
        O.write("sample_id,total,mapped,% mapped\n")
        for result in results:
            O.write(f"{result['sample_id']},{result['total']},{result['mapped']},{result['% mapped']:.2f}\n")

def job_collate_fasta_consensus(
    input_directories: list[FullPath],
    output_fasta: FullPath
):
    with open(output_fasta, "w") as O:
        for dir in input_directories:
            fasta_path = os.path.join(dir, "final_consensus.fasta")
            sample_name = os.path.basename(dir)
            with pysam.FastaFile(fasta_path) as fasta:
                for record in fasta.references:
                    seq = fasta.fetch(record)
                    O.write(f">{sample_name} {record}\n{seq}\n")



def job_median_depth_samtools(
    input_bam: FullPath,
    output_json: FullPath
):
    """Calculate median read depth across the genome using samtools depth.
    
    Args:
        input_bam (FullPath): Path to input BAM file with reads mapped to reference.
        output_json (FullPath): Path where output JSON file with median depth will be written.
    Returns:
        None
    """
    logging.info(f"Calculating median read depth with samtools depth. Input BAM: {input_bam}, Output JSON: {output_json}")
    tmp_depth_file = f"{output_json}.depth"
    cmd = f"samtools depth {input_bam} > {tmp_depth_file}"
    run_cmd(cmd)
    depths = []
    for l in open(tmp_depth_file):
        row = l.strip().split("\t")
        depths.append(int(row[2]))
    median_depth = np.median(depths) if depths else 0
    with open(output_json, "w") as O:
        json.dump({"median_depth": median_depth}, O)


def job_read_length_histogram(
    input_fastq: FullPath,
    output_tsv: FullPath
):
    """Generate a read length histogram from a FASTQ file and save as TSV.
    
    Args:
        input_fastq (FullPath): Path to input FASTQ file with reads.
        output_tsv (FullPath): Path where output TSV file with read length histogram will be written. Columns: read_length, count.
    Returns:
        None
    """
    # bins of 50 bp up to 20 kb
    bins = list(range(0, 20001, 50))
    logging.info(f"Generating read length histogram from {input_fastq}. Output TSV: {output_tsv}")
    length_counts = {b: 0 for b in bins}
    with pysam.FastxFile(input_fastq) as fastq:
        for entry in fastq:
            length = len(entry.sequence)
            # find the nearest bin
            bin_length = min(bins, key=lambda x: abs(x - length))
            length_counts[bin_length] += 1

    with open(output_tsv, "w") as O:
        O.write("read_length\tcount\n")
        for length, count in sorted(length_counts.items()):
            length = str(length) if length < 20000 else "20000+"
            O.write(f"{length}\t{count}\n")


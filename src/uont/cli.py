"""
CLI-related functions for the uont pipeline
"""

import argparse
from copy import deepcopy
import csv
import logging
import os
import sys
import rich_argparse
from rich.logging import RichHandler
from types import SimpleNamespace
from typing import get_args, get_origin, Literal
import yaml
from .types import FullPath
from .utils import JobStatus
from .utils import check_cli_dependencies, DEFAULT_CLI_DEPENDENCIES, get_dorado_model_dir, get_plassembler_db_dir
from .utils import g

from .process import (
    process_collate_barcode_fastqs,
)
from .workflow import (
    wf_assemble,
    wf_collate_amplicon_results,
    wf_amplicon,
    run_configured_workflow,
    wf_scrub,

)
import uont

from uont import __version__ as uont_version

g['uont_version'] = uont_version

def _fullpath_constructor(loader: yaml.Loader, node: yaml.Node) -> FullPath:
    value = loader.construct_scalar(node)
    return file_path(value)

yaml.SafeLoader.add_constructor("!FullPath", _fullpath_constructor)

def load_yaml_config(config_path: str) -> dict:
    """Load a YAML configuration file and return it as a dictionary.

    Args:
        config_path (str): Path to the YAML configuration file.
    Returns:
        dict: The loaded configuration as a dictionary.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file '{config_path}' does not exist")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f) or {}
    if not isinstance(config, dict):
        raise ValueError("Configuration file must contain a mapping at the top level")
    return config

def load_path_list_file(path_list_file: str) -> list[str]:
    """Load newline-delimited paths from a text file.

    Empty lines and lines starting with '#' are ignored.
    """
    if not os.path.exists(path_list_file):
        raise FileNotFoundError(f"Input directories file '{path_list_file}' does not exist")

    paths: list[str] = []
    with open(path_list_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            paths.append(file_path(line))

    if not paths:
        raise ValueError(
            f"Input directories file '{path_list_file}' does not contain any usable paths"
        )
    return paths

def file_path(path: str) -> str:
    """Convert a string path to an absolute path.

    Args:
        path (str): The input path.

    Returns:
        str: The absolute path.
    """
    # if not os.path.exists(path):
    #     logging.error(f"The file \"{path}\" does not exist. Exiting.")
    #     quit(1)
    return os.path.abspath(path)

def update_args_from_config(args: argparse.Namespace, config: dict) -> argparse.Namespace:
    updated_args = deepcopy(args)
    for key, value in config.items():
        normalized_key = key.replace("-", "_")
        if f"--{normalized_key}" in sys.argv:
            logging.info(
                "Command-line argument for %s takes precedence over config file value %s",
                normalized_key,
                value,
            )
            continue
        if hasattr(updated_args, normalized_key):
            setattr(updated_args, normalized_key, value)
        else:
            logging.warning(
                "Config file contains key '%s' which does not correspond to any command-line argument. Ignoring entry.",
                normalized_key,
            )
    return updated_args

def get_genome_sizes() -> dict:
    """Load genome size estimates from the included CSV file.

    Returns:
        dict: A mapping of organism names to their estimated genome sizes in base pairs.
    """
    genome_sizes_path = os.path.join(os.path.dirname(__file__), "data", "genome-sizes.csv")
    genome_sizes = {}
    for row in csv.DictReader(open(genome_sizes_path)):
        organism = row["organism"].strip()
        size = int(row["size"])
        genome_sizes[organism] = size
    return genome_sizes

def initialise_tools(
    args: argparse.Namespace,
) -> SimpleNamespace:
    """Initialize tools configuration from command-line arguments.

    Args:
        args (argparse.Namespace): Parsed command-line arguments from argparse.

    Returns:
        SimpleNamespace: Container of tool configuration attributes.
    """
    tools = SimpleNamespace()
    for key, value in vars(args).items():
        if key.endswith("_tool"):
            setattr(tools, key.replace("_tool", ""), value)
    return tools




def configure_logging(debug: bool = False, log_file: str | None = None) -> None:
    """Configure console logging and optional file logging for the CLI."""
    level = logging.DEBUG if debug else logging.INFO
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    console_handler = RichHandler(rich_tracebacks=True)
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    root_logger.addHandler(console_handler)

    if log_file:
        log_file = os.path.abspath(log_file)
        log_dir = os.path.dirname(log_file)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
        )
        root_logger.addHandler(file_handler)

        logging.info("Writing logs to file: %s", log_file)

def cli_uont():
    """Main entry point for the uont CLI."""
    parent_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )


    parent_parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )
    parent_parser.add_argument(
        "--preserve-temp-dirs",
        action="store_true",
        help="Preserve temporary directories for debugging",
    )
    parent_parser.add_argument(
        "--log-file",
        type=file_path,
        help="Optional path to a log file. If provided, logs are written to both console and file.",
    )
    parent_parser.add_argument(
        "--config",
        type=str,
        help="Path to YAML configuration file with tool parameters. Command-line arguments will override config file settings.",
    )
    parent_parser.add_argument(
        "--test",
        type=str,
        help="Make a fake output file for workflow testing. Provide the path to the output file to create.",
    )
    parent_parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + uont_version,
        help="Show the version of uont",
    )
    parser = argparse.ArgumentParser(
        description="A pipeline to process ONT data to a polished assembly.",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    workflow_parser = subparsers.add_parser(
        "workflow",
        aliases=["wf"],
        help="Run workflow commands",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    workflow_subparsers = workflow_parser.add_subparsers(dest="workflow_command", help="Available workflows")


    ##############################
    ##### Workflow: prepare ######
    ##############################

    prepare_wf_parser = workflow_subparsers.add_parser(
        "prepare",
        help="Prepare and collate fastq files from barcoded samples",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    prepare_wf_parser.add_argument(
        "--source-dir",
        type=str,
        required=True,
        help="The directory containing the barcoded fastq files",
    )
    prepare_wf_parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="The directory to output the processed fastq files",
    )
    prepare_wf_parser.add_argument(
        "--sample-sheet",
        type=str,
        required=True,
        help="The sample sheet containing the lab_id and barcode information",
    )
    prepare_wf_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="If set, the commands to collate fastq files will be printed but not executed",
    )
    prepare_wf_parser.add_argument(
        "--force",
        action="store_true",
        help="If set, existing output files will be overwritten during fastq collation",
    )

    ##############################
    ##### Workflow: scrub ######
    ##############################

    scrub_wf_parser = workflow_subparsers.add_parser(
        "scrub",
        help="Run read scrubbing workflow to remove adapters and dehumanise reads",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    scrub_wf_parser.add_argument(
        "--input-reads",
        type=file_path,
        required=True,
        help="Input HTS (bam or fastq) file containing the raw reads",
    )
    scrub_wf_parser.add_argument(
        "--output-reads",
        type=file_path,
        required=True,
        help="Output HTS (bam or fastq) file for the scrubbed reads",
    )
    scrub_wf_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="The number of threads to use for read scrubbing",
    )
    scrub_wf_parser.add_argument(
        "--dehumanise",
        action="store_true",
        help="If set, the dehumanisation step will be performed",
    )
    scrub_wf_parser.add_argument(
        "--sequencing-kit",
        type=str,
        choices=["SQK-RBK114-96"],
        help="The sequencing kit used for the run (required if using dorado for adapter removal)",
    )
    scrub_wf_parser.add_argument(
        "--adapter-removal-tool",
        type=str,
        default="dorado",
        choices=["porechop","dorado"],
        help="The tool to use for adapter removal in the scrubbing workflow",
    )
    scrub_wf_parser.add_argument(
        "--write-fastq",
        action="store_true",
        help="If input is bam and this is set, the output reads will also be written in fastq format"
    )


    from uont.utils import get_available_assemblers
    ##############################
    ##### Workflow: assemble #####
    ##############################

    assemble_wf_parser = workflow_subparsers.add_parser(
        "assemble",
        help="Assemble trimmed reads into contigs",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    assemble_wf_parser.add_argument(
        "--input-reads", 
        type=file_path, 
        required=True, 
        help="Input reads file (bam or fastq)"
    )
    assemble_wf_parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for assembly results",
    )
    assemble_wf_parser.add_argument(
        "--adapter-removal-tool",
        type=str,
        default="dorado",
        choices=["porechop","dorado"],
        help="The tool to use for adapter removal",
    )
    assemble_wf_parser.add_argument(
        "--fastq-filter-tool",
        type=str,
        default="chopper",
        choices=["chopper"],
        help="The tool to use for fastq filtering",
    )
    assemble_wf_parser.add_argument(
        "--genome-size-estimation-tool",
        type=str,
        default="lrge",
        choices=["lrge","autocycler"],
        help="The tool to use for genome size estimation",
    )
    assemble_wf_parser.add_argument(
        "--read-downsampling-tool",
        type=str,
        default="filtlong",
        choices=["filtlong"],
        help="The tool to use for read downsampling",
    )
    assemble_wf_parser.add_argument(
        "--assemblers",
        type=str,
        nargs="+",
        default=("flye", "miniasm","nextdenovo", "raven","plassembler","myloasm"),
        choices=get_available_assemblers(),
        help="The assembler(s) to use for assembly",
    )
    assemble_wf_parser.add_argument(
        "--polishing-tool",
        type=str,
        default="dorado",
        choices=["medaka","dorado"],
        help="The tool to use for polishing the assembly",
    )
    assemble_wf_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="The number of threads to use for assembly",
    )
    assemble_wf_parser.add_argument(
        "--threads-per-assembly",
        type=int,
        default=1,
        help="The number of threads to use per assembly job (used in autocycler assembly)",
    )
    assemble_wf_parser.add_argument(
            "--parallel-assembly-jobs",
            type=int,
            default=4,
            help="The number of parallel assembly jobs to run (used in autocycler assembly)",
        )
    assemble_wf_parser.add_argument(
        "--assembly-timeout-seconds",
        type=float,
        default=None,
        help="Per-assembly timeout in seconds for each parallel autocycler helper job",
    )
    assemble_wf_parser.add_argument(
        "--min-read-depth",
        type=int,
        default=25,
        help="Minimum read depth for subsampling (used in autocycler assembly)",
    )
    assemble_wf_parser.add_argument(
        "--max-contigs",
        type=int,
        default=80,
        help="Maximum number of contigs allowed (used in autocycler assembly)",
    )
    assemble_wf_parser.add_argument(
        "--max-samples",
        type=int,
        default=4,
        help="Maximum number of samples to process in parallel (used in autocycler assembly)",
    )
    assemble_wf_parser.add_argument(
        "--min-read-length",
        type=int,
        default=1000,
        help="Minimum read length for filtering (used in fastq filtering step)",
    )
    assemble_wf_parser.add_argument(
        "--min-q-score",
        type=int,
        default=10,
        help="Minimum average read quality score for filtering (used in fastq filtering step)",
    )
    assemble_wf_parser.add_argument(
        "--bam-for-dorado",
        type=str,
        # subpress help
        help = argparse.SUPPRESS,
    )
    assemble_wf_parser.add_argument(
        "--organism",
        type=str,
        choices=get_genome_sizes().keys(),
        help="Organism name for genome size (if not provided, will be estimated from the data)",
    )
    assemble_wf_parser.add_argument(
        "--lab-id",
        type=str,
        help="Lab ID to use for naming the final assembly and downsampled reads",
    )
    assemble_wf_parser.add_argument(
        "--link-id",
        type=str,
        help="Link ID to use for creating symbolic links to the final assembly",
    )
    assemble_wf_parser.add_argument(
        "--link-directory",
        type=str,
        help="Directory to create symbolic links for the final assembly",
    )
    assemble_wf_parser.add_argument(
        "--medaka-batch-size",
        type=int,
        default=100,
        help="Batch size for medaka polishing",
    )
    assemble_wf_parser.add_argument(
        "--models-directory",
        type=file_path,
        default=get_dorado_model_dir(),
        help="Directory containing custom Dorado models",
    )
    assemble_wf_parser.add_argument(
        "--save-filtered-reads",
        action="store_true",
        help="If set, the filtered reads from the assembly workflow will be saved to the output directory",
    )
    assemble_wf_parser.add_argument(
        "--save-unpolished-contigs",
        action="store_true",
        help="If set, the unpolished assembly will be saved to the output directory in addition to the polished assembly",
    )
    assemble_wf_parser.add_argument(
        "--additional-assemblies-dir",
        type=file_path,
        default=None,
        help="Optional directory containing additional assemblies to include in the final assembly set",
    )


    ########## END Workflow: assemble ##########

    


    ########### END Workflow subparsers ###########


    
    # Set up process subparser

    process_parser = subparsers.add_parser(
        "process",
        help="Run individual processes",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    process_subparsers = process_parser.add_subparsers(dest="process_command", help="Available processes")

    deps_parser = subparsers.add_parser(
        "check-dependencies",
        aliases=["deps"],
        help="Check whether required CLI tools are available on PATH",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    deps_parser.add_argument(
        "--tools",
        nargs="+",
        default=None,
        help=(
            "Optional list of tools to check. If omitted, checks defaults: "
            + ", ".join(DEFAULT_CLI_DEPENDENCIES)
        ),
    )

    setup_parser = subparsers.add_parser(
        "setup",
        help="Setup required tools and databases for uont",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )


    # Set up job and process subparsers

    job_parser = subparsers.add_parser(
        "job",
        help="Run individual jobs",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    job_subparsers = job_parser.add_subparsers(dest="job_command", help="Available jobs")
    

    # Set up custom parsers

    ####################################
    ##### Job: autocycler-assemble #####
    ####################################

    job_subparser = job_subparsers.add_parser(
        "autocycler-assemble",
        help="Run the complete autocycler assembly workflow",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    job_subparser.add_argument(
        "--input-fastq", 
        type=str, 
        required=True, 
        help="Input fastq file"
    )
    job_subparser.add_argument(
        "--output-dir", 
        type=str, 
        required=True, 
        help="Output directory for assembly results"
    )
    job_subparser.add_argument(
        "--genome-size", 
        type=int, 
        required=True, 
        help="Estimated genome size in base pairs"
    )
    job_subparser.add_argument(
        "--threads", 
        type=int, 
        default=4, 
        help="The number of threads to use"
    )
    job_subparser.add_argument(
        "--assembler", nargs="+", 
        type=str, 
        default=["miniasm"], 
        choices=["miniasm", "flye", "raven"], 
        help="Assembler to use"
    )
    job_subparser.add_argument(
        "--min-read-depth", 
        type=int, 
        default=25, 
        help="Minimum read depth for subsampling"
    )
    job_subparser.add_argument(
        "--max-contigs", 
        type=int, 
        default=80, 
        help="Maximum number of contigs allowed"
    )



    import inspect
    from . import process as process_module
    process_functions = [func for func in dir(process_module) if func.startswith("process_")]

    for func_name in process_functions:
        func = getattr(process_module, func_name)
        sig = inspect.signature(func)
        func_name = func_name.replace("_", "-").replace("process-", "")

        process_subparser = process_subparsers.add_parser(
            func_name,
            help=f"Run the {func_name} process",
            parents=[parent_parser],
            formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
        )

        for param in sig.parameters.values():

            if param.name=='kwargs':
                continue

            arg_name = f"--{param.name.replace('_', '-')}"
            arg_type = param.annotation if param.annotation is not param.empty else str

            arg = {
                'dest': param.name,
                'help': f"{param.name}",
            }
            if param.default is param.empty:
                arg['required'] = True
            else:
                arg['default'] = param.default
            if get_origin(arg_type) is Literal:
                arg['type'] = str
                arg['choices'] = get_args(arg_type)
            elif arg_type == bool:
                arg['action'] = 'store_true' if param.default is False else 'store_false'
                # remove the default value for boolean flags since it's implied by the action
                if 'default' in arg:                    
                    del arg['default']
            elif isinstance(arg_type, type):
                arg['type'] = arg_type
            elif arg_type == FullPath:
                arg['type'] = file_path
            else:
                arg['type'] = str

            process_subparser.add_argument(arg_name,**arg)


    # Add more job and process subparsers here following the same pattern 

    ## Get the names of all the job_subparsers added
    job_subparser_names = [name for name in job_subparsers.choices.keys()]
    ## Find all jobs in jobs.py 

    ## Find all functions in jobs.py that start with "job_" and add them as subcommands to the CLI
    from . import jobs as jobs_module
    
    job_functions = [func for func in dir(jobs_module) if func.startswith("job_")]

    for func_name in job_functions:
        func = getattr(jobs_module, func_name)
        sig = inspect.signature(func)
        if func_name.replace("job_", "").replace("_", "-") in job_subparser_names:
            continue

        job_subparser = job_subparsers.add_parser(
            func_name.replace("job_", "").replace("_", "-"),
            help=f"Run the {func_name} job",
            parents=[parent_parser],
            formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
        )
        for param in sig.parameters.values():

            if param.name=='kwargs':
                continue

            arg_name = f"--{param.name.replace('_', '-')}"
            arg_type = param.annotation if param.annotation is not param.empty else str
            arg = {
                'dest': param.name,
                'help': f"{param.name}",
            }

            if param.default is param.empty:
                arg['required'] = True
            else:
                arg['default'] = param.default
            if get_origin(arg_type) is Literal:
                arg['type'] = str
                arg['choices'] = get_args(arg_type)
            elif arg_type == list[FullPath]:
                arg['type'] = file_path
                arg['nargs'] = '+'
            elif arg_type == bool:
                arg['action'] = 'store_true' if param.default is False else 'store_false'
                # remove the default value for boolean flags since it's implied by the action
                if 'default' in arg:                    
                    del arg['default']
            elif isinstance(arg_type, type):
                arg['type'] = arg_type
            elif arg_type == FullPath:
                arg['type'] = file_path
            else:
                arg['type'] = str
            job_subparser.add_argument(arg_name,**arg)

    
    args = parser.parse_args()

    g['input_command'] = " ".join(sys.argv)

    configure_logging(debug=args.debug, log_file=args.log_file)
    if args.preserve_temp_dirs:
        uont.PRESERVE_TEMP_DIRS = True

    config_data = load_yaml_config(args.config) if args.config else None

    if config_data:
        args = update_args_from_config(args, config_data)

    if args.command in ("workflow", "wf"):
        if args.workflow_command == "prepare":
            process_collate_barcode_fastqs(
                args.source_dir, 
                args.output_dir, 
                args.sample_sheet,
                args.dry_run,
                args.force
            )

        elif args.workflow_command == "scrub":
            tools = initialise_tools(args)
            wf_scrub(
                input_reads=args.input_reads,
                output_reads=args.output_reads,
                tools=tools,
                threads=args.threads,
                dehumanise=args.dehumanise,
                sequencing_kit=args.sequencing_kit,
                write_fastq=args.write_fastq
            )

        elif args.workflow_command == "assemble":
            tools = initialise_tools(args)

            if args.organism:
                genome_size = get_genome_sizes().get(args.organism)
                logging.debug(f"Using genome size of {genome_size} for organism '{args.organism}'")
            else:
                genome_size = None
            if args.link_id and not args.link_directory:
                logging.error("If --link-id is provided, --link-directory must also be provided. Exiting.")
                sys.exit(1)
            wf_assemble(
                input_reads=args.input_reads,
                output_dir=args.output_dir,
                tools=tools,
                threads=args.threads,
                threads_per_assembly=args.threads_per_assembly,
                assemblers=args.assemblers,
                parallel_assembly_jobs=args.parallel_assembly_jobs,
                assembly_timeout_seconds=args.assembly_timeout_seconds,
                min_read_depth=args.min_read_depth,
                max_contigs=args.max_contigs,
                max_samples = args.max_samples,
                min_read_length=args.min_read_length,
                min_q_score=args.min_q_score,
                genome_size=genome_size,
                lab_id=args.lab_id,
                link_id=args.link_id,
                link_directory=args.link_directory,
                batch_size=args.medaka_batch_size,
                models_directory=args.models_directory,
                save_filtered_reads=args.save_filtered_reads,
                save_unpolished_contigs=args.save_unpolished_contigs,
                additional_assemblies_dir=args.additional_assemblies_dir
            )
            
        
        elif args.workflow_command == "run-config-workflow":
            run_configured_workflow(config_data)
        else:
            workflow_parser.print_help()
    
    elif args.command == "job":

        if args.job_command:

            from . import jobs as jobs_module
            func_name = f"job_{args.job_command}".replace("-", "_")
            if hasattr(jobs_module, func_name):

                func = getattr(jobs_module, func_name)
                import inspect
                sig = inspect.signature(func)
                kwargs = {
                    param.name: getattr(args, param.name)
                    for param in sig.parameters.values()
                    if param.name != "kwargs"
                }

                func(**kwargs)
            else:
                job_parser.print_help()
        else:
            job_parser.print_help()
    
    elif args.command == "process":
        if args.process_command:
            from . import process as process_module
            func_name = f"process_{args.process_command}".replace("-", "_")

            if hasattr(process_module, func_name):
                func = getattr(process_module, func_name)
                import inspect
                sig = inspect.signature(func)
                kwargs = {
                    param.name: getattr(args, param.name)
                    for param in sig.parameters.values()
                    if param.name != "kwargs"
                }
                func(**kwargs)
                
            else:
                process_parser.print_help()
        else:
            process_parser.print_help()
    elif args.command in ("check-dependencies", "deps"):
        available, missing = check_cli_dependencies(args.tools)

        if missing:
            logging.error("Missing dependencies: %s", ", ".join(missing))
            sys.exit(1)
        logging.info("All checked dependencies are available.")
    elif args.command == "setup":
        from .utils import setup_dorado, setup_plassembler_db
        setup_dorado()
        setup_plassembler_db()
    else:
        parser.print_help()

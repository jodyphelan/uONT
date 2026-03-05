"""
CLI-related functions for the uONT pipeline
"""

import argparse
import logging
import os
import rich_argparse
from types import SimpleNamespace
from typing import get_args, get_origin, Literal

from .process import (
    process_collate_barcode_fastqs,
)
from .workflow import (
    wf_assemble,
)

from .types import FullPath

def file_path(path: str) -> str:
    """Convert a string path to an absolute path.

    Args:
        path (str): The input path.

    Returns:
        str: The absolute path.
    """
    if not os.path.exists(path):
        logging.error(f"The file \"{path}\" does not exist. Exiting.")
        quit(1)
    return os.path.abspath(path)


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


def cli_uONT():
    """Main entry point for the uONT CLI."""
    parent_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )


    parent_parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )

    parser = argparse.ArgumentParser(
        description="A pipeline to process ONT data to a polished assembly.",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    workflow_parser = subparsers.add_parser(
        "workflow",
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

    ##############################
    ##### Workflow: assemble #####
    ##############################

    assemble_wf_parser = workflow_subparsers.add_parser(
        "assemble",
        help="Assemble trimmed reads into contigs",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    assemble_wf_parser.add_argument("--input", type=str, required=True, help="Input fastq file")
    assemble_wf_parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for assembly results",
    )
    assemble_wf_parser.add_argument(
        "--sample-name",
        type=str,
        help="The name of the sample being processed (used for naming output files)",
    )
    assemble_wf_parser.add_argument(
        "--adapter-removal-tool",
        type=str,
        default="porechop",
        choices=["porechop"],
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
        default="autocycler",
        choices=["autocycler"],
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
        "--assembler-tool",
        type=str,
        default="autocycler",
        choices=["autocycler","flye"],
        help="The assembler to use for assembly",
    )
    assemble_wf_parser.add_argument(
        "--polishing-tool",
        type=str,
        default="medaka",
        choices=["medaka"],
        help="The tool to use for polishing the assembly",
    )
    assemble_wf_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="The number of threads to use for assembly",
    )
    assemble_wf_parser.add_argument(
        "--min-read-depth",
        type=int,
        default=10,
        help="Minimum read depth for subsampling (used in autocycler assembly)",
    )
    assemble_wf_parser.add_argument(
        "--max-contigs",
        type=int,
        default=80,
        help="Maximum number of contigs allowed (used in autocycler assembly)",
    )

    ########### END Workflow subparsers ###########


    
    # Set up process subparser

    process_parser = subparsers.add_parser(
        "process",
        help="Run individual processes",
        parents=[parent_parser],
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
    )
    process_subparsers = process_parser.add_subparsers(dest="process_command", help="Available processes")




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
        "--sample-name", 
        type=str, 
        required=True, 
        help="The name of the sample being processed"
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
        default=10, 
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
            elif isinstance(arg_type, type):
                arg['type'] = arg_type
            elif arg_type == FullPath:
                arg['type'] = file_path
            else:
                arg['type'] = str
            job_subparser.add_argument(arg_name,**arg)

    
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.command == "workflow":
        if args.workflow_command == "prepare":
            process_collate_barcode_fastqs(
                args.source_dir, 
                args.output_dir, 
                args.sample_sheet
            )
        elif args.workflow_command == "assemble":
            tools = initialise_tools(args)
            wf_assemble(
                args.input,
                args.output_dir,
                tools,
                args.threads,
                args.min_read_depth,
                args.max_contigs,
            )
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
                result = func(**kwargs)
                if result is not None:
                    print(result)
            else:
                process_parser.print_help()
        else:
            process_parser.print_help()
    else:
        parser.print_help()

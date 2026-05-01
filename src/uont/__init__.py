"""
A package to process ONT data
to a polished assembly
"""

__version__ = "0.7.3"

import logging
from rich.logging import RichHandler

# Configure logging with rich
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[RichHandler(rich_tracebacks=True)]
)

# Import from jobs module
from .jobs import (
    job_fastq_filter_chopper,
    job_remove_adapters_porechop,
    job_downsample_filtlong,
    job_assemble_autocycler,
    job_polish_medaka,
)

# Import from process module
from .process import (
    Sample,
    process_collate_barcode_fastqs,
    process_fastq_filter,
    process_remove_adapters,
    process_estimate_genome_size,
    process_polish,
)

# Import from workflow module
from .workflow import (
    wf_assemble,
    make_dir_if_not_exists,
    run_configured_workflow,
)

# Import from cli module
from .cli import (
    cli_uONT,
    # cli_prepare,
    initialise_tools,
)

__all__ = [
    "__version__",
    "Sample",
    "process_collate_barcode_fastqs",
    "job_fastq_filter_chopper",
    "process_fastq_filter",
    "job_downsample_filtlong",
    "job_remove_adapters_porechop",
    "process_remove_adapters",
    "process_estimate_genome_size",
    "job_assemble_autocycler",
    "job_polish_medaka",
    "process_polish",
    "cli_uONT",
    "cli_prepare",
    "wf_assemble",
    "initialise_tools",
    "make_dir_if_not_exists",
    "run_configured_workflow",
]

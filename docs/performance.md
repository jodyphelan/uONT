# Performance

In this section, we provide some benchmarks for the computational performance of the uONT assembly workflow. We ran the workflow on a dataset of ONT reads from a bacterial genome, using an Apple MacbookPro M5. The `assemble` workflow was run using the the following command:

```bash
uONT workflow assemble --input ~/test-data/fastq/$1.fastq.gz --output $1 --threads 4 | tee $1.log
```



## Computational performace

### Overall runtime

This plot shows the total runtime of the uONT assembly workflow for multiple runs on the different samples. The total runtime is measured from the start of the workflow to the completion of the polishing step. This gives an overall sense of how long it takes to run the complete assembly process using uONT.

{{ plot_total_time() }}

### Task runtime breakdown

This plot shows the runtime of each individual task in the uONT assembly workflow, averaged across multiple runs. The tasks are ordered by the order of execution in the workflow. This breakdown helps identify which steps in the workflow are the most computationally intensive and may benefit from optimization.

{{ plot_execution_times() }}

### Comparison to old pipeline

This plot compares the total runtime of the uONT assembly workflow to a previous version of the pipeline that used different tools and processes. The comparison highlights any improvements in runtime performance achieved with the new workflow.

{{ plot_pipeline_runtimes() }}

This plot shows the same information but normalized to the runtime of the old pipeline, giving a clear visualization of the speedup achieved with the new workflow.

{{ plot_pipeline_speedup() }}

### Task runtime vs input file size

This plot shows the relationship between the input file size and the runtime of each task in the uONT assembly workflow. It helps to understand how the workflow scales with increasing data size and identify potential bottlenecks.

{{ plot_time_by_size() }}

## Biological performance

### Length of the final assembly

{{ plot_assembly_length() }}

### Number of contigs in the final assembly

{{ plot_num_contigs() }}

### N50 of final assembly

{{ plot_n50() }}

### GC content of final assembly

{{ plot_gc_content() }}
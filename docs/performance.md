# Computational performance

In this section, we provide some benchmarks for the computational performance of the uONT assembly workflow. We ran the workflow on a dataset of ONT reads from a bacterial genome, using an Apple MacbookPro M5. The `assemble` workflow was run using the the following command:

```bash
uONT workflow assemble --input ~/test-data/fastq/$1.fastq.gz --output $1 --threads 4 | tee $1.log
```



# Performance results

## Task runtime breakdown

This plot shows the runtime of each individual task in the uONT assembly workflow, averaged across multiple runs. The tasks are ordered by the order of execution in the workflow. This breakdown helps identify which steps in the workflow are the most computationally intensive and may benefit from optimization.

{{ plot_execution_times() }}


## Task runtime vs input file size

This plot shows the relationship between the input file size and the runtime of each task in the uONT assembly workflow. It helps to understand how the workflow scales with increasing data size and identify potential bottlenecks.

{{ plot_time_by_size() }}
"""
Basic example of a Mkdocs-macros module
"""

import math
import yaml
import os
import pandas as pd
import altair as alt

@alt.theme.register("large_font", enable=True)
def _large_font_theme() -> alt.theme.ThemeConfig:
    return alt.theme.ThemeConfig({
        "config": {
            "axis":   {"labelFontSize": 14, "titleFontSize": 16},
            "title":  {"fontSize": 18},
            "legend": {"labelFontSize": 14, "titleFontSize": 16},
            "header": {"labelFontSize": 14, "titleFontSize": 16},
        }
    })




def define_env(env):
    """
    This is the hook for defining variables, macros and filters

    - variables: the dictionary that contains the environment variables
    - macro: a decorator function, to declare a macro.
    - filter: a function with one of more arguments,
        used to perform a transformation
    """
    @env.macro
    def plot_execution_times():
        df = pd.read_csv('docs/assets/execution_times.csv')
        # Keep a stable, explicit ordering of categories on the x-axis.
        task_order = sorted(df["task"].dropna().unique().tolist())
        chart = (
            alt.Chart(df)
            .mark_boxplot()
            .encode(
                x=alt.X("task:N", sort=task_order, title="Task"),
                y=alt.Y("time_percent:Q", title="Time spent on task (%)"),
                color=alt.Color("task:N", legend=None),
                tooltip=["task:N"],
            )
            .properties(width=700, height=420, title="Time Spent per Task (%)")
            .interactive()
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"

    @env.macro
    def plot_time_by_size():
        df = pd.read_csv('docs/assets/execution_times.csv')
        df = df.sort_values(["file_size", "task"])
        # plot line of best fit for each task
        colours = {'01_job_remove_adapters_porechop': '#7fc97f',
        '02_job_fastq_filter_chopper': '#beaed4',
        '03_job_estimate_genome_size_lrge': '#fdc086',
        '04_job_downsample_filtlong': '#ffff99',
        '05_job_assemble_autocycler': '#386cb0',
        '06_job_polish_medaka': '#f0027f'}

        chart = (
            alt.Chart(df)
            .mark_point()
            .encode(
                x=alt.X("file_size:Q", title="Input file size (MB)"),
                y=alt.Y("seconds:Q", title="Time spent on task (seconds)"),
                color=alt.Color("task:N", title="Task", scale=alt.Scale(domain=list(colours.keys()), range=list(colours.values()))),
                tooltip=["task:N", "file_size:Q", "seconds:Q"],
            )
            .properties(width=700, height=420, title="Task Runtime vs Input File Size with Line of Best Fit")
            .interactive()
        )
        for task in df["task"].unique():
            task_df = df[df["task"] == task]
            line_of_best_fit = alt.Chart(task_df).mark_line(color=colours[task]).encode(
                x=alt.X("file_size:Q"),
                y=alt.Y("seconds:Q"),
            ).transform_regression("file_size", "seconds", method="linear")
            chart += line_of_best_fit

        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"
    

    @env.macro
    def plot_total_time():
        df = pd.read_csv('docs/assets/execution_times.csv')
        file_sizes = df.groupby("id")["file_size"].first().reset_index()
        file_sizes
        # calculate total time per log file
        total_times = df.groupby("id")["seconds"].sum().reset_index()
        total_times
        newdf = pd.merge(file_sizes, total_times, on="id")
        # get minutes
        newdf["minutes"] = newdf["seconds"] / 60
        # altair plot
        import altair as alt




        chart = (alt.Chart(newdf)
            .mark_point()
            .encode(
                x=alt.X("file_size:Q", title="File Size"),
                y=alt.Y("minutes:Q", title="Time (minutes)"),
                tooltip=["id:N", "file_size:Q", "minutes:Q"]
            )
            .properties(width=700, height=420, title="Execution Time vs File Size")
            .interactive()
            )
        
        line_of_best_fit = alt.Chart(newdf).mark_line(color="red").encode(
            x=alt.X("file_size:Q"),
            y=alt.Y("minutes:Q"),
        ).transform_regression("file_size", "minutes", method="linear")
        chart += line_of_best_fit
        chart_json = chart.to_json()
        return f"```vegalite\n{chart_json}\n```"

    @env.macro
    def plot_pipeline_runtimes():
        df = pd.read_csv('docs/assets/total_execution_times.csv')
        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("id:N", title="Run ID"),
                y=alt.Y("seconds:Q", title="Total Execution Time (seconds)"),
                color=alt.Color("pipeline:N", title="Pipeline"),
                xOffset=alt.XOffset("pipeline:N"),
                tooltip=["id:N", "pipeline:N", "seconds:Q"],
            )
            .properties(width=700, height=420, title="Total Execution Time by Pipeline")
            .interactive()
        )
        chart_json = chart.to_json()
        return f"```vegalite\n{chart_json}\n```"

    @env.macro
    def plot_pipeline_speedup():

        df = pd.read_csv('docs/assets/total_execution_times.csv')
        speedup_df = (
            df
            .pivot(index="id", columns="pipeline", values="seconds")
            .dropna(subset=["old-nanopore-wf", "uONT"])
            .reset_index()
        )

        speedup_df["speedup"] = speedup_df["old-nanopore-wf"] / speedup_df["uONT"]
        speedup_df["speedup_label"] = speedup_df["speedup"].map(lambda value: f"{value:.1f}x")

        chart = (
            alt.Chart(speedup_df)
            .mark_bar()
            .encode(
                x=alt.X("id:N", title="Run ID"),
                y=alt.Y("speedup:Q", title="Speedup vs old pipeline (x)"),
                tooltip=[
                    alt.Tooltip("id:N", title="Run ID"),
                    alt.Tooltip("old-nanopore-wf:Q", title="old-nanopore-wf (s)", format=".2f"),
                    alt.Tooltip("uONT:Q", title="uONT (s)", format=".2f"),
                    alt.Tooltip("speedup:Q", title="Speedup", format=".2f"),
                ],
            )
            .properties(width=700, height=420, title="Speedup vs old pipeline")
            .interactive()
        )
        
        
        chart_json = chart.to_json()
        return f"```vegalite\n{chart_json}\n```"

    @env.macro
    def plot_n50():
        df = pd.read_csv('docs/assets/assembly_qc.csv')

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("SampleID:N", title="Sample ID"),
                y=alt.Y("N50:Q", title="N50"),
                color=alt.Color("pipeline:N", title="Pipeline"),
                xOffset=alt.XOffset("pipeline:N"),
                tooltip=["SampleID:N", "pipeline:N", "N50:Q"],
            )
            .properties(width=700, height=420, title="Number of contigs by pipeline")
            .interactive()
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"
    
    @env.macro
    def plot_num_contigs():
        df = pd.read_csv('docs/assets/assembly_qc.csv')

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("SampleID:N", title="Sample ID"),
                y=alt.Y("Num_Contigs:Q", title="Number of contigs"),
                color=alt.Color("pipeline:N", title="Pipeline"),
                xOffset=alt.XOffset("pipeline:N"),
                tooltip=["SampleID:N", "pipeline:N", "Num_Contigs:Q"],
            )
            .properties(width=700, height=420, title="Number of contigs by pipeline")
            .interactive()
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"

    @env.macro
    def plot_gc_content():
        df = pd.read_csv('docs/assets/assembly_qc.csv')

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("SampleID:N", title="Sample ID"),
                y=alt.Y("GC_Content:Q", title="Number of contigs"),
                color=alt.Color("pipeline:N", title="Pipeline"),
                xOffset=alt.XOffset("pipeline:N"),
                tooltip=["SampleID:N", "pipeline:N", "GC_Content:Q"],
            )
            .properties(width=700, height=420, title="Number of contigs by pipeline")
            .interactive()
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"
    
    @env.macro
    def plot_assembly_length():
        df = pd.read_csv('docs/assets/assembly_qc.csv')

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("SampleID:N", title="Sample ID"),
                y=alt.Y("Length:Q", title="Assembly Length (bp)"),
                color=alt.Color("pipeline:N", title="Pipeline"),
                xOffset=alt.XOffset("pipeline:N"),
                tooltip=["SampleID:N", "pipeline:N", "Length:Q"],
            )
            .properties(width=700, height=420, title="Assembly length by pipeline")
            .interactive()
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"
    

    
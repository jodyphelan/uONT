"""
Basic example of a Mkdocs-macros module
"""

import math
import yaml
import os
import pandas as pd
import plotly.express as px
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
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"

    @env.macro
    def plot_time_by_size():
        df = pd.read_csv('docs/assets/execution_times.csv')
        df = df.sort_values(["file_size", "task"])
        chart = (
            alt.Chart(df)
            .mark_line(point=True)
            .encode(
                x=alt.X("file_size:Q", title="Input file size (MB)"),
                y=alt.Y("seconds:Q", title="Time spent on task (seconds)"),
                color=alt.Color("task:N", title="Task"),
                tooltip=["task:N", "file_size:Q", "seconds:Q"],
            )
            .properties(width=700, height=420, title="Task Runtime vs Input File Size")
        )
        json_text = chart.to_json()
        return f"```vegalite\n{json_text}\n```"
    

    
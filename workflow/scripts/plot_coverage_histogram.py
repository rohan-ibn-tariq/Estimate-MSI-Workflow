#!/usr/bin/env python3
"""
Plot histogram of MS region coverage fractions.
"""


import polars as pl
import altair as alt
import sys


input_file = snakemake.input.coverage
output_file = snakemake.output.plot
title = snakemake.params.title

df = pl.read_csv(
    input_file,
    separator="\t",
    has_header=False,
    new_columns=["chr", "start", "end", "name", "overlap_count", 
                 "bases_covered", "length", "fraction"],
    schema_overrides={
        "chr": pl.String,
        "start": pl.Int64,
        "end": pl.Int64,
        "name": pl.String,
        "overlap_count": pl.Int64,
        "bases_covered": pl.Int64,
        "length": pl.Int64,
        "fraction": pl.Float64,
    }
)

chart = alt.Chart(df).mark_bar().encode(
    x=alt.X('fraction:Q', bin=True, title='Coverage Fraction'),
    y=alt.Y('count()', title='Number of MS Regions')
).properties(
    width=700,
    height=400,
    title=title
)

chart.save(output_file)

print(f"Histogram saved to: {output_file}", file=sys.stderr)

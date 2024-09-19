from snakemake.utils import Paramspace
import pandas as pd

# declare a dataframe to be a paramspace
paramspace = Paramspace(pd.read_csv("params.csv", index_col='Unnamed: 0'))


rule all:
    input:
        # Aggregate over entire parameter space (or a subset thereof if needed)
        # of course, something like this can happen anywhere in the workflow (not
        # only at the end).
        expand("results/{params}.pickle", params=paramspace.instance_patterns)

rule simulate:
    output:
        # format a wildcard pattern like "alpha~{alpha}/beta~{beta}/gamma~{gamma}"
        # into a file path, with alpha, beta, gamma being the columns of the data frame
        f"results/{paramspace.wildcard_pattern}.pickle"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"alpha": ..., "beta": ..., "gamma": ...})
        simulation=paramspace.instance
    script:
        "simulate.py"

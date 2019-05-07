# Utilities for data pre-processing and plotting
This folder provides Python3 scripts for common data processing operations and plotting.

## Preprocessing
* `bin_sample_times.py` : reduces the number of time points by binning.


## Plotting
The script `plot_Q.py` plots stacked bar plots for ancestry estimates while (optionally) mantaining colors across K. A example is provided in the script `plot_demo.sh`

`plot_Q.py` has three required arguments:

1. the filepath to the Q matrix to plot  
2. population labels of each sample (one per line, corresponding to each row of Q)
3. and the order to plot populations (one per line, corresponding to each row of Q)

Calling

```
python plot_Q.py Q3 samplelabels poporder
```

will plot the results for 3 populations, where samples are grouped by their populations defined in `samplelabels` and ordered according to `poporder.` Within each population, samples are sorted by the ancestry proportions of their major contributing population. Calling
```
python plot_Q.py -h
```
prints a full list of plotting options to the console. A subset of options are explained in detail below.


### Matching Colors
To match labels across runs, `plot_Q.py` provies the optional argument `--match-Q`. `--match-Q` takes the filepath of a folder containing multiple `Q` matrices each assumed to be named `QK` for each value of `K`. By specifying `--match-Q`, `plot_Q.py` internally re-arranges the order of populations sequentially for each value of `K` to match colors across runs.

For example, suppose we ran `K=2,3,4` and placed the Q matrices `Q2, Q3, Q4` in the folder `./Q_matrices.` We now want to plot the results for each `K`, matching the colors across runs. This can be achieved as follows

```
python plot_Q.py ./Q_matrices/Q2 samplelabels poporder
python plot_Q.py ./Q_matrices/Q3 samplelabels poporder --match_Q ./Q_matrices 
python plot_Q.py ./Q_matrices/Q4 samplelabels poporder --match_Q ./Q_matrices 

```
You can alternatively match colors across different runs by specifying a folder of Q matrices from a previous run.


### Plotting Subsets
`plot_Q.py` also supports plotting subsets of samples using the optional argument `--subset`. `--subset` gives the path to a file containing a list of populations (one per line) in the order to be plot.

Example:
```
python plot_Q.py ./Q_matrices/Q3 samplelabels poporder --subset subset
```


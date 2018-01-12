# Introduction
Dystruct is a model-based approach for inferring ancestry proportions from time-series genotype data, in particular ancient DNA. Dystruct corrects for genetic drift between ancient and contemporary populations

## Getting Started

### Quick Start

```
git clone https://github.com/tyjo/dystruct
cd dystruct
make
./bin/dystruct
./run.sh
```

### OpenMP
Dystruct depends on the OpenMP API for parallel computation. Unfortunately, OpenMP is not supported on the clang compiler shipped with macOS. Mac users will need to compile Dystruct using a different C++ compiler. To install gcc using homebrew, simply call

```
brew install gcc6
```

If your compiler does not support OpenMP, you will likely see the error

```
clang: error: unsupported option '-fopenmp'
```

You can either remove this flag from the makefile, which disables parallel computation, or install a different compiler per the above directions.


## Running Dystruct

### Parallel Computation
Dystruct parallelizes certain sections of its algorithm for faster computation using the OpenMP library. To specify the number of threads OpenMP looks for the environment variable OMP\_NUM\_THREADS. Thus, if you want to use 2 threads, you should call

```
export OMP_NUM_THREADS=2
```

before running.



### Arguments
Dystruct takes several command line arguments:

```
Usage:   dystruct [options]

Options:
	--input FILE                Genotype file path. An LOCI x INDIVIDUAL matrix of genotypes. The header
                                is the sample time in generations. Samples must be ordered in increasing
                                generation time.
	--output STR                A file prefix for output files.
	--npops INT                 Number of populations.
	--loci INT                  Number of loci. This should match the number of loci in the input file.
	--pop-size INT              Effective population size for all populations.
	--seed INT                  Random seed used to initialize variational parameters
	--hold-out-fraction DOUBLE  (=0) Optional. Partitions nloci * hold_out_fraction loci into a hold out
                                set. The hold out set contains at most one site per individual.
	--hold-out-seed INT         (=28149) Optional. Random seed used to partition SNP data into hold out
                                and training sets. Use the same seed across replicates to fix the hold
                                out set.
	--step-size-power DOUBLE    (=-0.6) Optional. Adjusts step size for stochastic variational inference.
                                step_size = (iteration - offset)^step_power after the first 10000
                                iterations. The offset ensures the step size does jump between iteration
                                10000 and 10001. Must be in [-1,-0.5).
	--labels FILE               Optional. Experimental. Population label file path for supervised analysis.
                                Labels should be in {0,...,npops - 1}. One label per line in the same order
                                as the input matrix. Individuals without a population assignment should be
                                labeled by -1.
```

For convenience, a shell script that runs Dystruct on a test data set and sets the appropriate environment variables is provided (run.sh).

#### Input Files
Dystruct takes as input a whitespace delimited file of genotypes. Each column is a vector of genotypes (0, 1, or 2) for a single sampled individual, where the header is the time (in generations) when the individual was sampled relative to the first sample. Sites are assumed to be biallelic, where a 0 denotes that an individual a homozygote with respect to one of the alleles, a 1 denotes that an individual is a heterozygote, and a 2 denotes that an individual is homozygous for the alternate allele. Missing sites should be labled by 9. **The columns and are assumed to be in ascending order of sample time.** Each row is a locus. Below is an example for 2 sampled loci from 3 individuals at time 0 and 1.

```
0 0 1
0 1 2
0 1 2
2 2 0
```
There are two individuals sampled a time 0: the first with genotypes (0, 0, 2) and the second with genotypes (1, 1, 2).

An example file is provided in

```
./supp/data/BASELINE_100GEN_10000LOCI_120D_01/samples
```


#### Output Files
Dystruct outputs two files with point estimates for inferred parameters, and a temporary file to monitor convergence. The prefix of these files is specified through the --output command line argument. The inferred parameter files are:

- freqs : the inferred allele frequencies at all time steps
- theta : the inferred ancestry proportions for all samples

Dystruct also outputs a temporary file, temp\_theta, with the current estimates of the variational parameters for ancestry proportions. Each row is an individual, and each column a population. The number in each entry refers to the number of loci currently assigned to a population. Dystruct uses a relatively strick convergence criterion, and this file can be used to moniter looser convergence criteria.

#### Hold-out Fraction and Choosing Parameters
Dystruct takes an option argument to hold out a subset of loci. After convergence, the program outputs the log likelihood of the held-out set, treating each time point independently. Specifically, the likelihood each binomial observation (genotype) is evaluated using the inferred allele frequencies at that time point and the individual's ancestry proportions. This has the benefit of being independent of the choice of effective population size, and can also be used to evaluate different choices of number of populations (*K*): choose the set of parameters with the highest held out log likelihood.

Nonetheless, in our experience choosing *K* is of minor importance. Instead, like related ancestry inference programs, Dystruct should be run over a range of *K* and each *K* analyzed separately. Each *K* can give insight to a different relationship between ancient and modern samples.

#### Labeled Data
Dystruct takes an optional argument with population assignment for ancient samples. These are treated as known and fixed. Ancestry is inferred for the remaining unlabeled samples. This feature is currently experimental, and has yet to be thoroughly investigated using simulations.

### Running Time
Running time per iteration depends on two factors: i) how many times the algorithm has previously seen a loci, and ii) number of time points. Earlier iterations tend to be much slower than later iterations because it takes longer for the allele frequency estimates to converge. The time spent per iteration decreases after the first two or three cycles through each loci.

The update step for allele frequency estimatation takes quadratic time with respect to the number of time points. Therefore it is desirable to group ancient samples into as few time points as possible. For example, ancient samples with overlapping confidence intervals for carbon date estimates should be combined into a single time point.

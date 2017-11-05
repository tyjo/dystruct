# Introduction
Dystruct is a model-based approach for inferring ancestry proportions from time-series genotype data, in particular ancient DNA.

# Compilation and Dependencies
Dystruct depends on several libraries in [Boost](http://www.boost.org), and OpenMP API for parallel computation. If Boost is not installed on your machine, you can download and install a local copy using the instructions here. Dystruct has been compiled and tested with Boost 1.62.0, however newer releases of Boost _should_ work as well. First, download a copy of the .tar.gz file to the dystruct directory from [https://sourceforge.net/projects/boost/files/boost/1.62.0/](https://sourceforge.net/projects/boost/files/boost/1.62.0/). Once the download has completed, move to this directory in the Terminal.

```
dystruct tyjo$ cd boost_1_62_0
boost_1_62_0 tyjo$ 
boost_1_62_0 tyjo$ ./bootstrap.sh --with-libraries=program_options --prefix=. --with-toolset=gcc
boost_1_62_0 tyjo$ ./b2 install
```


## Parallel computation
Dystruct parallelizes certain sections of its algorithm for faster computation using the OpenMP library. To specify the number of threads OpenMP looks for the environment variable OMP\_NUM\_THREADS. Thus, if you want to use 2 threads, you should call

```
export OMP_NUM_THREADS=2
```

before running.

## Mac users
OpenMP is not currently supported on the clang compiler that ships with macOS. In order to compile Dystruct, you will need to download a different C++ compiler. If you use homebrew, you can simply install g++ using

```
brew install gcc6
```
Then, after running the bootstrap.sh script above, edit the project-config.jam file. Replace the lines

```
 10 if ! gcc in [ feature.values <toolset> ]
 11 {
 12     using gcc ;
 13 }
```

with the lines

```
10 if ! gcc in [ feature.values <toolset> ]
11 {
12     using gcc : : g++-6 ; 
13 }
```

If you use the wrong compiler, you will most likely see the error

```
clang: error: unsupported option '-fopenmp'
```
## Compiling
Once the dependencies are satisfied, you can compile with

```
make
```

This outputs the Dystruct binary to bin. Run using

```
./bin/dystruct
```

If you get an error, you are most likely missing your Boost install from your LD\_LIBRARY\_PATH. If you installed Boost locally, you can update the environment variable:

```
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/PATH/TO/BOOST/boost_1_62_0/stage
```

# Running Dystruct
Dystruct takes several command line arguments:

```
Dystruct Program Usage:
  --help                              Print help message
  -i [ --input ] arg                  Genotype file path. An L x D matrix of 
                                      genotypes, where the header is the sample
                                      time in generations. Samples must be 
                                      ordered in increasing generation time.
  -o [ --output ] arg                 A file prefix that is appended to all 
                                      output files.
  -k [ --npops ] arg                  Number of populations.
  -l [ --nloci ] arg                  Number of loci.
  -z [ --pop_size ] arg               Specifies population size for all 
                                      populations.
  -s [ --seed ] arg                   Random seed used to initialize 
                                      variational parameters
  -f [ --hold_out_fraction ] arg (=0) Fraction of loci to hold out.
  -h [ --hold_out_seed ] arg (=28149) Random seed used to partition SNP data 
                                      into hold out and training sets. Use the 
                                      same seed across replicates to keep the 
                                      hold out set fixed.
  -b [ --labels ] arg                 Population labels for supervised 
                                      analysis. Labels should be in 
                                      {0,...,npops - 1}. One label per line in 
                                      the same order as the input matrix. 
                                      Individuals without a population 
                                      assignment should be labeled by -1.
```

If the hold\_out\_fraction is none zero, a fraction of all loci will be held out, and a lower bound on the posterior predictive distribution evaluated when the algorithm has converged. This gives an indication of how well Dystruct is performing on the data, and can be compared across runs to choose the best run for a particular value of K. Similarly, this can be compared for different values of K to choose the best K. The hold\_out\_seed provides a seed to the random number generator to select loci. Keep this value consistant, or leaving it set to the default, will ensure that the same loci are held out across runs.

For convenience, a shell script that runs Dystruct on a test data set and sets the appropriate environment variables is provided. A common complication is that your Boost installation is not in your LD\_LIBRARY\_PATH. If you get an error, you will need to set the appropriate LD\_LIBRARY\_PATH in this file.

## Input Files
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
supp/data/BASELINE\100GEN\_10000LOCI\_120D\_01
```


## Output Files
Dystruct outputs two files with point estimates for inferred parameters, and a temporary file to monitor convergence. The prefix of these files is specified through the --output command line argument. The inferred parameter files are:

- freqs : the inferred allele frequencies at all time stes
- theta : the inferred ancestry proportions for all samples

Dystruct also outputs a temporary file, temp\_theta, with the current estimates of the variational parameters for ancestry proportions. Each row is an individual, and each column a population. The number in each entry refers to the number of loci currently assigned to a population. Dystrut uses a relatively strick convergence criterion, and this file can be used to moniter looser convergence criteria.

## Monitering convergence
Dystruct outputs progress every 1000 iterations by displaying a global step size, and the average change in number of loci assigned to a population for each individual (appears as delta). Once the average loci assigned is less than 1 the algorithm terminates, as this is now beyond the precision available for point estimates.
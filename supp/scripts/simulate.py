#!/usr/bin/env/ python

# Simulate data under the Wright-Fisher model of genetic drift.
import argparse
import numpy as np
import pickle
import sys

parser = argparse.ArgumentParser(description='Simulate genotype data across generations.')
parser.add_argument('-k', type=int, help='number of populations')
parser.add_argument('-l', type=int, help='number of loci')
parser.add_argument('-s', '--seed', type=int, help='random seed to run simulation')
parser.add_argument('--set', type=str, help='simulation set')
parser.add_argument('--sizes', nargs='+', type=int,
                    help='effective population sizes for each population')
parser.add_argument('--samples', nargs='+', type=int,
                    help='number of individuals sampled at each generation')
args = parser.parse_args()
# global parameters
K = args.k
L = args.l
sim_set = args.set
seed = args.seed
pop_sizes = args.sizes
sample_times = args.samples
T = len(sample_times)
alpha = np.array([1.0 / K for i in range(K)])
outfile1 = "samples"
outfile2 = "theta"
outfile3 = "freqs"


def generate_allele_frequencies():
    beta = np.zeros((T, K, L))

    for k in range(K):
        for l in range(L):
            beta[0][k][l] = np.random.uniform(0.2, 0.8)
    
    for t in range(1, len(sample_times)):
        for k in range(K):
            for l in range(L):
                n = np.random.binomial(2*pop_sizes[k], beta[t-1][k][l])
                beta[t][k][l] = n/(2*pop_sizes[k])

    return beta


def generate_theta():
    theta = [ [ [ 0 for k in range(K) ] for j in range(sample_times[t]) ] for t in range(T) ]
    for t in range(T):
        for d in range(sample_times[t]):
            if sim_set.lower() == "baseline":
                dr = np.random.dirichlet(alpha).tolist()
                for k in range(K):
                    theta[t][d][k] = dr[k]
            elif sim_set.lower() == "clustered":
                if t == 0:
                    dr = np.random.dirichlet([20, 1, 1])
                elif t == len(sample_times)/2 - 1:
                    dr = np.random.dirichlet([1, 20, 1])
                elif t == len(sample_times) - 1:
                    dr = np.random.dirichlet([1, 1, 20])
                for k in range(K):
                    theta[t][d][k] = dr[k]
            elif sim_set.lower() == "continuous":
                for k in range(K):
                    theta[t][d][k] = 0.
                theta[t][d][np.random.randint(0,3)] = 1.
            # single ancestral population
            elif sim_set.lower() == "ancestral":
                if t == 0:
                    theta[t][d][0] = 1
                    theta[t][d][1] = 0
                elif t == len(sample_times) - 1:
                    dr = np.random.dirichlet([1,1])
                    for k in range(K):
                        theta[t][d][k] = dr[k]
            else:
                print("bad simulation set")
                sys.exit(1)
    return theta


def generate_samples(theta, beta):
    X = [ [ [ 0 for l in range(L) ] for j in range(sample_times[t]) ] for t in range(T) ]

    for t in range(T):
        for d in range(sample_times[t]):
            for l in range(L):
                p = 0.0
                for k in range(K):
                    p += theta[t][d][k]*beta[t][k][l]
                if p >= 1 and abs(1 - p) < 0.0001:
                    p = 1
                X[t][d][l] = np.random.binomial(2, p)
    return X

if __name__ == "__main__":
    print('Simulating:')
    print('   {} individuals'.format(sum(sample_times)))
    print('   {} loci'.format(L))
    print('   {} generations'.format(len(sample_times)))
    print('   {} populations'.format(K))
    #print('   {} population sizes'.format(pop_sizes))
    #print('   {} sample times'.format(sample_times))
    np.random.seed(seed)

    beta = generate_allele_frequencies()
    theta = generate_theta()
    X = generate_samples(theta, beta)

    SNP = [] # generation sampled + SNPs
    theta_fmt = [] # generation sampled + theta
    for t in range(T):
        for d in range(len(X[t])):
            
            row = [t]
            for l in X[t][d]:
                row.append(l)
            SNP.append(row)

            row = []
            for k in theta[t][d]:
                row.append(k)
            theta_fmt.append(row)

    SNP = np.array(SNP, dtype=int)
    np.savetxt(outfile1, SNP.T, fmt="%d", delimiter="\t")

    theta_fmt = np.array(theta_fmt)
    np.savetxt(outfile2, theta_fmt, fmt="%.5f", delimiter="\t")

    with open(outfile3, "w") as f:
        for t in range(beta.shape[0]):
            if (sample_times[t] > 0):
                f.write(str(t) + "\n");
                for l in range(beta.shape[2]):
                    for k in range(beta.shape[1] - 1):
                        f.write(str(beta[t][k][l]) + "\t")
                    f.write(str(beta[t][K-1][l]) + "\n")
                f.write("\n")
        

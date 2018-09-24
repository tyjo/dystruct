import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Bin generation times to reduce number of time points. Does not change the range of sample times.")
parser.add_argument("input", help="Path to file with original generation times. Output file is saved to INPUT-reduced.")
parser.add_argument("--bucket_size", default=50, help="Width of each bucket. Divides generation times into TOTAL_TIME/BUCKET_SIZE buckets\n" + \
                                                    "Where the n-th bucket, starting from 0, contains individuals within the range\n" + \
                                                    "[n*BUCKET_SIZE, (n+1)*BUCKET_SIZE). Generation times are replace by the mean\n" + \
                                                    "within each bucket.")
args = parser.parse_args()
input_file = args.input

time_gen = open(input_file, "r").readlines()
time_gen = [float(gen.strip("\n")) for gen in time_gen]
time_gen = np.array(time_gen)

# make sure generations begin at 0
time_gen = time_gen - time_gen.min()

bucket_size = 50
buckets = [[] for i in range(int(time_gen.max() / bucket_size)+1)]
for idx,g in enumerate(time_gen):
    bucket_idx = int(g / bucket_size) 
    
    if g == time_gen.max():
        buckets[-1].append(g)
    elif g == 0:
        buckets[0].append(g)
    else:
        buckets[bucket_idx].append(g)

bucket_idx = 0
for idx,g in enumerate(time_gen):
    bucket_idx = int(g / bucket_size)

    if g == time_gen.max() or g == 0:
        continue
    else:
        time_gen[idx] = np.mean(buckets[bucket_idx])

np.savetxt(input_file + "-reduced", time_gen, fmt="%i")

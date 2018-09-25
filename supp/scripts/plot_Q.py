"""
Copyright (C) 2017-2018 Tyler Joseph <tjoseph@cs.columbia.edu>

This file is part of Dystruct.

Dystruct is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Dystruct is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Dystruct.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import scipy.cluster.hierarchy
import scipy.cluster.vq
import scipy.stats
from scipy import spatial
from scipy.linalg import norm
from itertools import permutations
from bokeh.charts import Bar, output_file, show

def match_one_Q(QK, Qprv):
    """Match the columns (populations) of Q to Qprv.
    """
    npops = QK.shape[1]
    distances = np.zeros((npops, npops))
    for id1, row1 in enumerate(QK.T):
        for id2, row2 in enumerate(Qprv.T):
            distances[id1, id2] = np.square(row1 - row2).sum()
    distances[:,npops-1] = 1E10
    assignments = [i for i in range(npops)]
    for i in range(npops):
        label = np.unravel_index(np.nanargmin(distances), (npops, npops))
        distances[label[0]] = np.array([np.nan for i in range(npops)])
        distances[:,label[1]] = np.array([np.nan for i in range(npops)])
        assignments[label[1]] = label[0]
    return assignments


def match_Q(Q, K, matchQ_path):
    """Match the populations for current Q to
    previous runs.
    """

    # find the lowest contiguous value of K to
    # begin matching colors
    min_k = K
    for k in range(K, 1, -1):
        try:
            Qmin = np.genfromtxt("Q" + str(k))
            min_k = k
        except:
            break
    print("matching colors from K =", min_k, "to", K)

    # sequentially match columns for each
    # value of Q
    Qprv = np.genfromtxt(matchQ_path + "/Q" + str(min_k))
    for i in range(min_k, K+1):
        Qnxt = np.genfromtxt(matchQ_path + "/Q" + str(i))
        order = match_one_Q(Qnxt, Qprv)
        Qnxt = Qnxt[:,order]
        Qprv = np.copy(Qnxt)
    order = match_one_Q(Q, Qprv)
    Q = Q[:,order]
    return Q


def plot_k(Q, labels, order, fontsize, spacing):
    z = np.zeros(Q.shape[1])
    groups = [ [z for i in range(spacing)] for pop in order ]
    labels_grouped = [ ["" for i in range(spacing)] for pop in order]
    for label,mixt in zip(labels, Q):
        groups[order.index(label)].append(mixt)
        labels_grouped[order.index(label)].append(label)

    for idx,group in enumerate(groups):
        group_z = group[:7]
        group_r = group[7:]
        if len(group_r) <= 1:
            continue

        # group individuals by major population contributor
        group_r = np.array(group_r)
        ids = np.argmax(group_r, axis=1)
        group_r = group_r[np.argsort(ids)]

        # sort within each population cluster
        ids = np.argmax(group_r, axis=1)
        sort_by = scipy.stats.mode(ids)[0][0]
        for idy in ids:
            rows = group_r[idy == ids].tolist()
            rows.sort(key = lambda row: row[sort_by])
            group_r[idy == ids] = np.array(rows)
        groups[idx] = group_z + group_r.tolist()

    labels_ordered = []
    prev_label = "-1"
    for group in labels_grouped:
        for label in group:
            if label == prev_label:
                label = ""
            else:
                prev_label = label
            labels_ordered.append(label)
    Q = np.zeros((Q.shape[0] + 7*len(groups), Q.shape[1]))
    idx = 0
    for g in groups:
        for m in g:
            Q[idx] = m
            idx += 1

    ancestry_matrix = Q

    plt.figure()
    fig, ax = plt.subplots(nrows=1, ncols=1)

    m = matplotlib.cm.get_cmap("tab20")
    colors = [m(9*i % 20) for i in range(20)]
    ancestry_matrix = pd.DataFrame(ancestry_matrix)
    ancestry_matrix.plot.bar(ax=ax, stacked=True, legend=False, width=1, linewidth=0.01, figsize=(25,5), color=colors)#colormap=matplotlib.cm.get_cmap("tab20"))

    ax.set_yticklabels([])
    ax.set_xticklabels(labels_ordered, fontsize=fontsize)#, rotation=45)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.patch.set_visible(False)
    fig.patch.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout()

    k = Q.shape[1]
    plt.savefig("dystruct_k" + str(k) + ".png", transparent=True, dpi=1000, pad_inches=0)
    plt.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Visualize ancestry estimates.")
    parser.add_argument("filepath",   type=str,                help="Path to Q matrix to plot.")
    parser.add_argument("K",          type=int,                help="Value of K to plot.")
    parser.add_argument("samplelabels",  type=str,             help="File path to population labels for each sample (one per line).")
    parser.add_argument("poporder",   type=str,                help="File path giving the order to plot each population (one per line).")
    parser.add_argument("--match-Q",  type=str,  default=None, help="Match colors across K. The argument to --match-Q is a path to a folder" +
                                                                   " of Q matrices, each named Qk for each value of K.")
    parser.add_argument("--subset",   type=str, default=None,  help="List of a subset of populations (one per line) from pop order to plot.")
    parser.add_argument("--fontsize", type=int, default=6,     help="Size of population labels.")
    parser.add_argument("--spacing",  type=int, default=7,     help="Size of spaces between populations.")
    args = parser.parse_args()

    path = args.filepath
    K = args.K
    samplelabels = args.samplelabels
    poporder = args.poporder
    matchQ = args.match_Q
    subset = args.subset
    fontsize = args.fontsize
    spacing = args.spacing

    Q = np.genfromtxt(path)

    assert Q.shape[1] == K, "error: value of K does not match Q matrix"

    if matchQ is not None:
        Q = match_Q(Q, K, matchQ)

    print("plotting...")

    samplelabels = open(samplelabels, "r").readlines()
    samplelabels = [label.strip("\n") for label in samplelabels]
    poporder = open(poporder, "r").readlines()
    poporder = [pop.strip("\n") for pop in poporder]

    if subset is not None:
        print("using subset of individuals in", subset)
        subset = open(subset, "r").readlines()
        subset = [pop.strip("\n") for pop in subset]

        sublabels = []
        subids = []
        for idx, label in enumerate(samplelabels):
            if label in subset:
                sublabels.append(label)
                subids.append(idx)

        samplelabels = sublabels
        poporder = subset
        Q = Q[subids]

    plot_k(Q, samplelabels, poporder, fontsize, spacing)




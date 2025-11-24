#!/usr/bin/env python3
"""
Plot sensor model, target positions, estimates and confidence ellipses.

Usage:
  python3 tools/plot_estimates.py --sensors sensors.csv --targets targets.csv \
      --estimates estimates.csv --cov covariances.npy

Inputs:
  - `sensors.csv` rows: id,x,y
  - `targets.csv` rows: id,x,y
  - `estimates.csv` rows: id,x,y
  - `covariances.npy`: numpy array of shape (N,2,2) or a single 2x2 covariance

The script draws 95% confidence ellipses for each covariance.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import csv
import math


def read_xy_csv(path):
    rows = []
    with open(path, newline='') as f:
        r = csv.reader(f)
        for line in r:
            if not line:
                continue
            # accept id,x,y or x,y
            if len(line) >= 3:
                _, x, y = line[:3]
            else:
                x, y = line[:2]
            rows.append((float(x), float(y)))
    return np.array(rows)


def confidence_ellipse(cov, mean, n_std=2.4477, **kwargs):
    # n_std = 2.4477 approximates 95% for 2D (sqrt(chi2.ppf(0.95,2)))
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    width, height = 2 * n_std * np.sqrt(vals)
    angle = math.degrees(math.atan2(vecs[1, 0], vecs[0, 0]))
    return Ellipse(xy=mean, width=width, height=height, angle=angle, **kwargs)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sensors', required=True)
    p.add_argument('--targets', required=False)
    p.add_argument('--estimates', required=False)
    p.add_argument('--cov', required=False)
    p.add_argument('--show', action='store_true')
    p.add_argument('--outfile', default='artifacts/plot.png')
    args = p.parse_args()

    fig, ax = plt.subplots(figsize=(8, 8))

    if args.sensors:
        s = read_xy_csv(args.sensors)
        ax.scatter(s[:, 0], s[:, 1], marker='^', color='C2', label='Sensors')
        for i, (x, y) in enumerate(s):
            ax.annotate(f'S{i}', (x, y), textcoords='offset points', xytext=(4, 4))

    if args.targets:
        t = read_xy_csv(args.targets)
        ax.scatter(t[:, 0], t[:, 1], marker='x', color='C1', label='Targets')

    covs = None
    if args.cov:
        covs = np.load(args.cov)
        # normalize shape
        if covs.ndim == 2 and covs.shape == (2, 2):
            covs = np.expand_dims(covs, 0)

    if args.estimates:
        e = read_xy_csv(args.estimates)
        ax.scatter(e[:, 0], e[:, 1], marker='o', color='C0', label='Estimates')
        for i, (x, y) in enumerate(e):
            ax.annotate(f'E{i}', (x, y), textcoords='offset points', xytext=(4, -8))
        if covs is not None:
            # if counts mismatch, either reuse first cov or clip
            n = len(e)
            if covs.shape[0] < n:
                covs = np.vstack([covs] + [covs[0]] * (n - covs.shape[0]))
            for i in range(n):
                cov = covs[i]
                mean = (e[i, 0], e[i, 1])
                ell = confidence_ellipse(cov, mean, edgecolor='C0', facecolor='none', lw=1.5)
                ax.add_patch(ell)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend()
    ax.grid(True)
    # ensure output directory exists
    import os
    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.outfile, dpi=200)
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()

#!/usr/bin/env python

"""plot_raster.py

Contour plot of a raster.

"""
import os
import argparse
import glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# parse args
parser = argparse.ArgumentParser(description="Contour plot of a raster")
parser.add_argument('filename', help="File containing Raster to plot")
parser.add_argument('-l', dest='logscale', action='store_true', help='Plot on logarithmic scale')
parser.add_argument('-s', dest='save', action='store_true', help="Save the image to file")
parser.add_argument('-o', dest='offline', action='store_true', help="Don't display interactive plot")
args = parser.parse_args()
fn = args.filename
assert os.path.exists(fn), "File does not exist: {}".format(fn)

# load data, first read nodata value
with open(fn) as fh:
    for i in range(5):
        fh.readline()
    nodata = float(fh.readline().split()[-1])

# now read data
data = np.loadtxt(fn, skiprows=6)

# fix nodata
if np.any(data == nodata):
    data[data == nodata] = np.nan

# now plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)
if args.logscale:
    plt.contourf(data, cmap='seismic', locator=matplotlib.ticker.LogLocator())
else:
    plt.contourf(data, cmap='seismic')
cbar = plt.colorbar()
plt.gca().set_aspect('equal', adjustable='box')
plt.title(fn)
plt.tight_layout()
if args.save:
    plt.savefig(fn + ".png")
if not args.offline:
    plt.show()

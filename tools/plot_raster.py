#!/usr/bin/env python

"""plot_raster.py

Contour plot of a raster.

"""
import os
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def plot_file(fn, args):
    print("Plotting file: {}".format(fn))

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

    # make the plot
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

    # save to file?
    if args.save:
        savefn = fn + ".png"
        print("Saving plot to: {}".format(savefn))
        plt.savefig(savefn)

    # view interactive?
    if not args.offline:
        plt.show()

    # clear figure
    plt.clf()


def main():
    # parse args
    parser = argparse.ArgumentParser(description="Create contour plots of raster files (multiple files can be specified)")
    parser.add_argument('filename', nargs='+', help="File containing Raster to plot")
    parser.add_argument('-l', '--logscale', action='store_true', help='Plot on logarithmic scale (default is linear)')
    parser.add_argument('-s', '--save', action='store_true', help="Save the image to file (default is not to save to file)")
    parser.add_argument('-o', '--offline', action='store_true', help="Don't display interactive plot (default is to display the interactive plot)")
    args = parser.parse_args()

    # plot all files specified
    for fn in args.filename:
        plot_file(fn, args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""plot_raster.py

Contour plot of a raster.

"""
import os
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def load_raster(fn):
    with open(fn) as fh:
        for i in range(5):
            fh.readline()
        nodata = float(fh.readline().split()[-1])

    # now read data
    data = np.loadtxt(fn, skiprows=6)

    # fix nodata
    if np.any(data == nodata):
        data[data == nodata] = np.nan

    return data


def plot_file(fn, args, data_diff):
    print("Plotting file: {}".format(fn))

    # load data, first read nodata value
    data = load_raster(fn)

    # are we plotting differences
    if data_diff is not None:
        data -= data_diff

    # plot absolute values
    if args.abs:
        data = np.absolute(data)

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
        savefn = os.path.basename(fn) + ".png"
        print("Saving plot to: {}".format(savefn))
        plt.savefig(savefn)

    # view interactive?
    if not args.offline:
        plt.show()

    if args.offline:
        # clear figure
        plt.clf()
        plt.close()


def main():
    # parse args
    parser = argparse.ArgumentParser(description="Create contour plots of raster files (multiple files can be specified)")
    parser.add_argument('filename', nargs='+', help="File containing Raster to plot")
    parser.add_argument('-l', '--logscale', action='store_true', help='Plot on logarithmic scale (default is linear)')
    parser.add_argument('-s', '--save', action='store_true', help="Save the image to file (default is not to save to file)")
    parser.add_argument('-o', '--offline', action='store_true', help="Don't display interactive plot (default is to display the interactive plot)")
    parser.add_argument('-d', '--diff', default=None, help="Plot the difference from the file specified by this option (default is not to plot differences)")
    parser.add_argument('-a', '--abs', action='store_true', help="Plot absolute values (default is not to take absolute values)")
    args = parser.parse_args()

    # are we plotting differences
    if args.diff is not None:
        print("Plotting differences of Rasters relative to '{}''".format(args.diff))
        data_diff = load_raster(args.diff)
    else:
        data_diff = None

    # absolute values
    if args.abs:
        print("Plotting absolute values")

    # plot all files specified
    for fn in args.filename:
        plot_file(fn, args, data_diff)


if __name__ == "__main__":
    main()

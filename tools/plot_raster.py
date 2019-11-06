#!/usr/bin/env python

"""plot_raster.py

Contour plot of a raster.

"""
import os
import argparse
import subprocess

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
        plt.contourf(data, cmap='seismic', locator=matplotlib.ticker.LogLocator(), vmin=args.min, vmax=args.max)
    else:
        plt.contourf(data, cmap='seismic', vmin=args.min, vmax=args.max)
    cbar = plt.colorbar()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(fn)
    plt.tight_layout()

    # save to file?
    savefn = None
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

    return savefn


def main():
    # parse args
    parser = argparse.ArgumentParser(description="Create contour plots of raster files (multiple files can be specified)")
    parser.add_argument('filename', nargs='+', help="File containing Raster to plot")
    parser.add_argument('-l', '--logscale', action='store_true', help='Plot on logarithmic scale (default is linear)')
    parser.add_argument('-s', '--save', action='store_true', help="Save the image to file (default is not to save to file)")
    parser.add_argument('-o', '--offline', action='store_true', help="Don't display interactive plot (default is to display the interactive plot)")
    parser.add_argument('-d', '--diff', default=None, help="Plot the difference from the file specified by this option (default is not to plot differences)")
    parser.add_argument('-a', '--abs', action='store_true', help="Plot absolute values (default is not to take absolute values)")
    parser.add_argument('-m', '--min', type=float, default=None, help="Min value to use in the colour map (default is the minimum value)")
    parser.add_argument('-u', '--max', type=float, default=None, help="Max value to use in the colour map (default is the maximum value)")
    parser.add_argument('-g', '--gif', action='store_true', help="Create gif from saved images (requires ImageMagick convert command line tool)")
    parser.add_argument('-D', '--delay', type=int, default=20, help="Argument to pass to `convert -delay` (default is 20)")
    args = parser.parse_args()

    if args.offline:
        matplotlib.use('Agg')

    # are we plotting differences
    if args.diff is not None:
        print("Plotting differences of Rasters relative to '{}'".format(args.diff))
        data_diff = load_raster(args.diff)
    else:
        data_diff = None

    # absolute values
    if args.abs:
        print("Plotting absolute values")

    # plot all files specified
    save_filenames = []
    for fn in args.filename:
        savefn = plot_file(fn, args, data_diff)
        save_filenames.append(savefn)

    # create gif animation
    if args.gif and args.save:
        prefix = os.path.commonprefix(save_filenames).split("_")[0]
        if not len(prefix):
            prefix = "output"
        gifname = prefix + ".gif"
        cmd = "convert -delay {} {} {}".format(args.delay, " ".join(save_filenames), gifname)
        print("Creating animation: {}".format(gifname))
        subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":
    main()

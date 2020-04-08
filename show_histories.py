#!/usr/bin/python3
# -*- coding: utf-8 -*-
__version__ = '0.0.7' # Time-stamp: <2020-02-10T19:21:06Z>

"""Show histories."""

import matplotlib.pyplot as plt
import argparse
import traceback


parser = argparse.ArgumentParser()
parser.add_argument("specs", metavar="LABEL=CSV", nargs='+', type=str)
parser.add_argument("--ax1", default="supply_of_labors", type=str)
parser.add_argument("--ax2", default=None, type=str)
parser.add_argument("--ylabel1", default=None, type=str)
parser.add_argument("--ylabel2", default=None, type=str)
parser.add_argument("--xlabel", default=None, type=str)
parser.add_argument("--title", default="Supply of Labors", type=str)

ARGS = parser.parse_args()


def save_history (path, history):
    import csv
    with open(path, 'w') as f:
        keys = list(history.keys())
        epochs = len(history[keys[0]])
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        r = [dict([[k, history[k][i]] for k in keys])
             for i in range(epochs)]
        writer.writerows(r)

def load_history (path):
    import csv
    with open(path) as f:
        reader = csv.DictReader(f)
        h = {}
        for k in reader.fieldnames:
            h[k] = []
        for row in reader:
            for k in reader.fieldnames:
                h[k].append(float(row[k]))
        return h


def run():
    specs = []
    epochs = []
    for s in ARGS.specs:
        label, fname = s.split('=', 1)
        history = load_history(fname)
        epochs.append(len(history[ARGS.ax1]))
        s = {'label': label, 'ax1': history[ARGS.ax1]}
        if ARGS.ax2 is not None:
            s['ax2'] = history[ARGS.ax2]
        specs.append(s)

    epochs = range(max(epochs))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    color = 0
    for s in specs:
        ax1.plot(epochs, s['ax1'], 'C' + str(color), label=s['label'])
        color += 1

    if ARGS.ax2 is not None:
        ax2 = ax1.twinx()
        for s in specs:
            ax2.plot(epochs, s['ax2'], 'C' + str(color), label=s['label'])
            color += 1

    h1, l1 = ax1.get_legend_handles_labels()
    if ARGS.ax2 is not None:
        h2, l2 = ax2.get_legend_handles_labels()
    else:
        h2, l2 = [], []
    ax1.legend(h1 + h2, l1 + l2)

    if ARGS.title is not None:
        plt.title(ARGS.title)
    if ARGS.xlabel is not None:
        ax1.set_xlabel(ARGS.xlabel)
    if ARGS.ylabel1 is not None:
        ax1.set_ylabel(ARGS.ylabel1)
    if ARGS.ylabel2 is not None:
        ax2.set_ylabel(ARGS.ylabel2)

    plt.show()

if __name__ == "__main__":
    run()

#!/usr/bin/python3

import os
import sys
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import csv
from csv import writer, reader


def output_csv(csvf, strain, ind, rep, dialect, out="outfile.csv"):
    with open(csvf) as f:
        f.readline()
        csvreader = reader(f, dialect)
        with open(out, "a+") as o:
            csv_writer = writer(o)
            for row in csvreader:
                csv_writer.writerow([row[ind], strain, rep, row[0]])
    return

def create_plot(csvf):
    df = pd.read_csv(csvf)
    ax = sns.lineplot(x=df['Time'], y=df['OD'], hue=df['Treatment'], ci="sd")
    plt.ylabel("OD$_{600}$")
    plt.xlabel("Time (minutes)")
    plt.show()


if __name__ == "__main__":
    # txt file specifying which columns belong to which strain
    # format: "strain\tcol1\tcol2[...]"
    keys = sys.argv[1]
    # csv containing actual data, first row needs to be names
    # first column needs to be time, in minutes or other integer based increments
    csvf = sys.argv[2]
    # name of output file, if none defaults to "output.csv"
    try:
        out = sys.argv[3]
    except IndexError:
        out = "outfile.csv"

    with open(csvf) as f:
        csvdialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)

    strains = {}
    with open(keys) as f:
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        csvreader = reader(f, dialect)
        for strain in csvreader:
            strains[strain[0]] = strain[1:]

    # get col headers
    with open(csvf) as f:
        csvreader = reader(f, csvdialect)
        for header in csvreader:
            headers = header
            break

    # create outfile and populate new headers
    with open(out, "w+") as o:
        csv_writer = writer(o)
        csv_writer.writerow(["OD", "Treatment", "Replicate", "Time"])
    for strain, cols in strains.items():
        rep = 1
        for col in cols:
            ind = headers.index(col)
            output_csv(csvf, strain, ind, rep, csvdialect)
            rep += 1
    create_plot(out)

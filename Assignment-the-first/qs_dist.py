#!/usr/bin/env python

"""This program parses through a FASTQ file and generates a per base distribution of quality scores."""

import matplotlib.pyplot as plt
import Bioinfo
import gzip
import argparse

def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="A program to generate a distribution graph of FASTQ quality scores")
    parser.add_argument("-f", "--file", help="Input filename", required=True)
    parser.add_argument("-p", "--plot", help="Plot filename", required=True)
    return parser.parse_args()

# store the command line args in variables
args = get_args()
file= args.file
plot = args.plot

read_length = 0
linecount = 0

# Determine the length of the sequence lines to be used for list initialization 
with gzip.open(file, "rt") as fr:
    fr.readline()
    seqline = fr.readline()
    read_length = len(seqline)

# initialize list to store ongoing sums
my_list = [0] * (read_length-1)

# parse the FASTQ file and sum the quality scores at each base pair
with gzip.open(file, "rt") as fr:
    for line in fr:
        line = line.strip()     
        if linecount % 4 == 3:
            for i, q in enumerate(line):
                my_list[i] += Bioinfo.convert_phred(q) 
        linecount+=1

print("# Base Pair" + "\t" + "Mean Quality Score")
for i, qv in enumerate(my_list):
    my_list[i] = (qv / (linecount/4))
    print(str(i) + "\t" + str(my_list[i]))

# Generate plots of the quality score distribution
x = range(len(my_list))
y = my_list
plt.plot(x, y)
plt.xlabel("Base Pair Position")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Scores per Base Pair Position")
plt.savefig(plot)
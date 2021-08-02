#!/usr/bin/env python

""" This program performs demultiplexing of RNA seq results by removal of read pairs that have been index swapped or have unknown indexes
    and the organization of all the read pairs with matched indices into specific buckets corresponding to their index. Given a pair of 
    read files and index files, will separate each read pair according to which category it falls into: matched index pair, swapped index 
    pair, unknown index pair. """

import argparse
import Bioinfo
import gzip

def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="A program to demultiplex sequencing files")
    #parser.add_argument("-f", "--files", help="Input FASTQ filenames (one read pair, one index pair)", nargs='4', required=True)
    parser.add_argument("-i", "--indexes", help="Index filename", required=True)
    #parser.add_argument("-o", "--ofile", help="Output filename", required=True)
    #parser.add_argument("-d", "--dir", help="Output directory", required=True)
    return parser.parse_args()

# store the command line args in variables
args = get_args()
#files= args.files
indexes = args.indexes
#ofile = args.ofile
#dir = args.dir

# Dictionary to store list of index seqs and other related index info (sample, group, treatment, index)
# Key = index sequence, Value = Dictionary identifying sample, group, treatment, and index values
index_dict = {}

# extract the index sequences and other related info from text file and store in index_dict
with open(indexes, "rt") as fr:
    fr.readline() 
    for line in fr:
        line = line.strip()
        cols = line.split('\t')
        index_info_dict = {}
        for col in cols:
            index_info_dict['sample'] = cols[0]
            index_info_dict['group'] = cols[1]
            index_info_dict['treatment'] = cols[2]
            index_info_dict['index'] = cols[3]
        index_dict[cols[4]] = index_info_dict

#print(index_dict)

# Dictionary to store indexes and their reverse complements // to avoid calling rev_comp for every index read while parsing index 2 fq file
# Keys = Index, Value = Reverse Complement of Index 
rc_index_dict = {}
for i in index_dict:
    rc_index_dict[i]= Bioinfo.rev_comp(i)

#print(index_dict)
#print(rc_index_dict)

###############
# for testing, remove later
#index_dict = {"GTCC": "i1_ut_GTCC_", "AGAA": "i2_ut_AGAA_"}
###############

# open all the output FASTQ files (48 matched files + 2 swapped + 2 unknown) to write to when parsing the input FASTQ files
swapped_R1 = gzip.open("unit_output/swapped_R1.fq.gz", "wt")
swapped_R2 = gzip.open("unit_output/swapped_R2.fq.gz", "wt")
unknown_R1 = gzip.open("unit_output/unknown_R1.fq.gz", "wt")
unknown_R2 = gzip.open("unit_output/unknown_R2.fq.gz", "wt")
# Dictionary to store all the output file handlers to be able to reference when parsing through the input FASTQ files
# Keys = name of file handle (same as file name), Values = the file handle
file_handlers_dict = {"swapped_R1.fq.gz":swapped_R1, "swapped_R2.fq.gz":swapped_R2, "unknown_R1.fq.gz":unknown_R1, "unknown_R2.fq.gz":unknown_R2}
for i in index_dict:
    r1= i + "_" + Bioinfo.rev_comp(i) + "_R1.fq.gz" 
    r2 = i + "_" + Bioinfo.rev_comp(i) + "_R2.fq.gz"
    file_handlers_dict[r1] = gzip.open("unit_output/" + r1, "wt") 
    file_handlers_dict[r2] = gzip.open("unit_output/" + r2, "wt")


# would this be easier ?? than just storing the 26 counters instead of redundant 52?
# Dictionary to store counts for all matched indexes, swapped indexes, unknown indexes, and low quality indexes
# Key = Bucket, Value = count
counter_dict = {}
for i in file_handlers_dict:
    counter_dict[i] = 0

#print(counter_dict)

# need to use gzip for the actual files
# open input files (pair of read FASTQs and pair of index FASTQs)
fq1 = gzip.open("../TEST-input_FASTQ/unit_test_R1.fq.gz", "rt")
fq2 = gzip.open("../TEST-input_FASTQ/unit_test_R2.fq.gz", "rt")
fq3 = gzip.open("../TEST-input_FASTQ/unit_test_R3.fq.gz", "rt")
fq4 = gzip.open("../TEST-input_FASTQ/unit_test_R4.fq.gz", "rt")

# Begin parsing the 4 input files one read at a time in parallel
while True:
    r1_lines = [0] * 4
    r2_lines = [0] * 4
    i1_lines = [0] * 4
    i2_lines = [0] * 4
    for i in range(4):
        r1_lines[i] = fq1.readline().strip()
        r2_lines[i] = fq2.readline().strip()
        i1_lines[i] = fq3.readline().strip()
        i2_lines[i] = fq4.readline().strip()
    
    if r1_lines[i] == '':
        break
    
    # print(r1_lines)
    # print(r2_lines)
    # print(i1_lines)
    # print(i2_lines)

    index1_seq = i1_lines[1]
    index2_seq = i2_lines[1]
    index1_qs = i1_lines[3]
    index2_qs = i2_lines[3]

    bucket = ""
    # check if both indexes meet the quality score cutoff
    if Bioinfo.meets_Qcutoff(index1_qs) and Bioinfo.meets_Qcutoff(index2_qs):
        # check if both indexes are valid // don't contain N's
        if index1_seq in index_dict and index2_seq in rc_index_dict:
            # check if the 2 indexes are reverse comp of each other
            if index1_seq == rev_comp(index2_seq):
                bucket = index1_seq
                # increment counter
            else:
                bucket = "swapped"
                # increment counter
        else:
            bucket = "unknown"
            # increment counter
    else:
        bucket = "unknown"
        # increment counter

    print(bucket)

# close the input files 
fq1.close()
fq2.close()
fq3.close()
fq4.close()

# close all the output FASTQ files 
for file in file_handlers_dict:
    file_handlers_dict[file].close()







###
# Making sure input files are zipped and FASTQ files? do a check?
###

"""
    For the 2 read FASTQ files and 2 index FASTQ files, will extract one read at at time (line by line) and append the read to its appropriate bucket file
        1. Extract and store the header [only up to the space // all 4 files should have the same header up to the space]
            [will append the pair of index sequences to the end of the header before outputting to files]
        2. Extract and store the read 1 and read 2 sequence lines as well as index 1 and index 2 sequence lines 
        3. Extract and store the third line of read 1 and read 2 
        4. Extract and store the quality score lines of read1, read 2, and index 1, index 2 

    Temp variables to include: header, R1seq, R2seq, I1seq, I2seq, line3, R1q, R2q, I1q, I2q

        At this point you should know which bucket your reads should go to and you have all the info stored for both reads
            To the pair of R1 and R2 files for the appropriate bucket (files should already be open):
                output the modified header to each file
                output the seq line to each file
                output the 3rd line to each file
                output the read quality score to each file 

    When you reach the end of files, close out all the open output files 
    Output the counts for read-pairs with matched indexes (one per index-pair), with index-hopping, and with unknwon indexes
"""
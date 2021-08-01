#!/usr/bin/env python

""" This program performs demultiplexing of RNA seq results by removal of read pairs that have been index swapped or have unknown indexes
    and the organization of all the read pairs with matched indices into specific buckets corresponding to their index. Given a pair of 
    read files and index files, will separate each read pair according to which category it falls into: matched index pair, swapped index 
    pair, unknown index pair. """

import argparse
import Bioinfo

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
# Key = string of index sequence, Value = other related index info 
index_dict = {}

# extract the index sequences and other related info from text file and store in index_dict
with open(indexes, "r") as fr:
    fr.readline()
    for line in fr:
        line = line.strip()
        cols = line.split('\t')
        indexinfo = ""
        for col in cols:
            indexinfo = indexinfo + col + "_" 
        seq = cols[4]
        index_dict[seq] = indexinfo


###############
# for testing, remove later
index_dict = {"GTCC": "i1_ut_GTCC_", "AGAA": "i2_ut_AGAA_"}
###############


# open all the output FASTQ files (48 matched files + 2 swapped + 2 unknown) to write to when parsing the input FASTQ files
swapped_R1 = open("output/swapped_R1.fq", "w")
swapped_R2 = open("output/swapped_R2.fq", "w")
unknown_R1 = open("output/unknown_R1.fq", "w")
unknown_R2 = open("output/unknown_R2.fq", "w")
# Dictionary to store all the output file handlers to be able to reference when parsing through the input FASTQ files
# Keys = name of file handle (same as file name), Values = the file handle
file_handlers = {"swapped_R1":swapped_R1, "swapped_R2":swapped_R2, "unknown_R1":unknown_R1, "unknown_R2":unknown_R2}
for i in index_dict:
    ########
    # make all the keys the same as the file names? including .fq part?
    ########
    r1= "output/" + i + "_" + Bioinfo.rev_comp(i) + "_R1.fq" 
    r2 = "output/" + i + "_" + Bioinfo.rev_comp(i) + "_R2.fq"
    file_handlers[i] = open(r1, "w") 
    file_handlers[Bioinfo.rev_comp(i)] = open(r2, "w")



# open input files (pair of read FASTQs and pair of index FASTQs)
fq1 = open("../TEST-input_FASTQ/unit_test_R1.fq", "r")
fq2 = open("../TEST-input_FASTQ/unit_test_R2.fq", "r")
fq3 = open("../TEST-input_FASTQ/unit_test_R3.fq", "r")
fq4 = open("../TEST-input_FASTQ/unit_test_R4.fq", "r")


# close the input files 
fq1.close()
fq2.close()
fq3.close()
fq4.close()

# close all the output FASTQ files 
for file in file_handlers:
    file_handlers[file].close()


"""
    Once your output files are open and ready to write to, will begin parsing through the input FASTQ files for read1, read2, and index1, index2
        Will go line by line to compare index 1 with index 2
        For each pair of indices we will determine which bucket it should be sent to and output the 4 lines for each pair of reads to either the R1 or R2 file in corresponding bucket
        Before going through the input files, make sure to initialize counter variables for each bucket

    For the 2 read FASTQ files and 2 index FASTQ files, will extract one read at at time (line by line) and append the read to its appropriate bucket file
        1. Extract and store the header [only up to the space // all 4 files should have the same header up to the space]
            [will append the pair of index sequences to the end of the header before outputting to files]
        2. Extract and store the read 1 and read 2 sequence lines as well as index 1 and index 2 sequence lines 
        3. Extract and store the third line of read 1 and read 2 
        4. Extract and store the quality score lines of read1, read 2, and index 1, index 2 

    Temp variables to include: header, R1seq, R2seq, I1seq, I2seq, line3, R1q, R2q, I1q, I2q

        Set temporary variable to store the bucket that the reads will be sent to after going through the following conditions
            If index1 and index2 RC are in the index dictionary // could be either matched or swapped 
                If index1 and index2 both meet the quality score cutoff:
                    if index 1 equals index2 RC:
                        send to matched index bucket and increment counter 
                    Else swapped indexes
                        send to swapped bucket and increment counter
                Else if either or both indexes don't meet cutoff: 
                    send to trash and increment counter      
            Else if either one or both indexes not in the dictionary // means not valid index:
                send to trash bin and increment counter
        At this point you should know which bucket your reads should go to and you have all the info stored for both reads
            To the pair of R1 and R2 files for the appropriate bucket (files should already be open):
                output the modified header to each file
                output the seq line to each file
                output the 3rd line to each file
                output the read quality score to each file 

    When you reach the end of files, close out all the open output files 
    Output the counts for read-pairs with matched indexes (one per index-pair), with index-hopping, and with unknwon indexes
"""
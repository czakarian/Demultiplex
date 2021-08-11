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
    parser.add_argument("-f", "--files", help="Input FASTQ filenames (one read pair, one index pair)", nargs=4, required=True)
    parser.add_argument("-i", "--indexes", help="Index filename", required=True)
    parser.add_argument("-d", "--direct", help="Output directory", required=True)
    parser.add_argument("-c", "--counts", help="Counts filename", required=True)
    return parser.parse_args()

# store the command line args in variables
args = get_args()
files= args.files
indexes = args.indexes
direct = args.direct
countsfile = args.counts

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

# Dictionary to store indexes and their reverse complements // to avoid calling rev_comp for every index read while parsing index 2 fq file
# Keys = Index, Value = Reverse Complement of Index 
rc_index_dict = {}
for i in index_dict:
    rc_index_dict[Bioinfo.rev_comp(i)]= i

# Dictionary to store counts for all matched indexes read-pairs
# Key = matched index pair, Value = count
counter_matched_dict = {}
for i in index_dict:
    counter_matched_dict[i] = 0

# Dictionary to store counts for all swapped index read-pairs
# Key = swapped index pair, value = count
counter_swapped_dict = {}

# variable to store number of unknown read-pairs
counter_unknown = 0
counter_swapped = 0
counter_matched = 0

# open all the output FASTQ files (48 matched files + 2 swapped + 2 unknown) to write to when parsing the input FASTQ files // open as zipped files
#make sure to have gzip, "wt" and .fq.gz
swapped_R1 = gzip.open(direct + "/swapped_R1.fq.gz", "wt")
swapped_R2 = gzip.open(direct + "/swapped_R2.fq.gz", "wt")
unknown_R1 = gzip.open(direct + "/unknown_R1.fq.gz", "wt")
unknown_R2 = gzip.open(direct + "/unknown_R2.fq.gz", "wt")

# Dictionary to store all the output file handlers to be able to reference when parsing through the input FASTQ files
# Keys = name of file handle (same as file name), Values = the file handle
file_handlers_dict = {"swapped_R1.fq.gz":swapped_R1, "swapped_R2.fq.gz":swapped_R2, "unknown_R1.fq.gz":unknown_R1, "unknown_R2.fq.gz":unknown_R2}
for i in index_dict:
    r1 = index_dict[i]["index"] + "_R1.fq.gz" 
    r2 = index_dict[i]["index"] + "_R2.fq.gz"
    file_handlers_dict[r1] = gzip.open(direct + "/" + r1, "wt") 
    file_handlers_dict[r2] = gzip.open(direct + "/" + r2, "wt")

# open input files (pair of read FASTQs and pair of index FASTQs)
fq1 = gzip.open(files[0], "rt")
fq2 = gzip.open(files[1], "rt")
fq3 = gzip.open(files[2], "rt")
fq4 = gzip.open(files[3], "rt")




# Begin parsing the 4 input files one read at a time in parallel
while True:
    # use a dict instead of list to store the 4 lines per each of the 4 reads??
    fq_reads = {"r1":{1:"", 2:"", 3:"", 4:""}, "r2":{1:"", 2:"", 3:"", 4:""}, "i1":{1:"", 2:"", 3:"", 4:""}, "i2":{1:"", 2:"", 3:"", 4:""}}
    for i in range(4):
        fq_reads["r1"][i] = fq1.readline().strip()
        fq_reads["r2"][i] = fq4.readline().strip()
        fq_reads["i1"][i] = fq2.readline().strip()
        fq_reads["i2"][i] = fq3.readline().strip()  
    
    # exit when you reach end of the file
    if fq_reads["r1"][0] == '':
        break

    # store the rc of index2 seq so you only need to call the actual fxn once
    rc = Bioinfo.rev_comp(fq_reads["i2"][1])

    # append i1 and i2 rc indexes to the header
    header = fq_reads["r1"][0].split(' ')[0] + "-" + fq_reads["i1"][1] + "-" + rc
    
    # determine which bucket your read will be sent to
    bucket = ""

    # check if both indexes are valid / don't contain N's and meet the quality score cutoff
    if fq_reads["i1"][1] in index_dict and fq_reads["i2"][1] in rc_index_dict and Bioinfo.meets_Qcutoff(fq_reads["i1"][3]) and Bioinfo.meets_Qcutoff(fq_reads["i2"][3]):
        # check if the 2 indexes are reverse comp of each other
        if fq_reads["i1"][1] == rc:
            bucket = index_dict[fq_reads["i1"][1]]["index"]
            counter_matched_dict[fq_reads["i1"][1]] += 1
            counter_matched += 1
        else:
            bucket = "swapped"
            counter_swapped += 1
            i1i2 = fq_reads["i1"][1] + "-" + rc
            i2i1 = rc + "-" + fq_reads["i1"][1]
            if i1i2 in counter_swapped_dict:
                counter_swapped_dict[i1i2] += 1
            elif i2i1 in counter_swapped_dict:
                counter_swapped_dict[i2i1] += 1
            else:
                counter_swapped_dict[i1i2] = 1          
    else:
        bucket = "unknown"
        counter_unknown += 1

    # append each read to the corresponding bucket file
    r1file = file_handlers_dict[bucket + "_R1.fq.gz"]
    r2file = file_handlers_dict[bucket + "_R2.fq.gz"]

    r1file.write(fq_reads["r1"][0] + "\n" + fq_reads["r1"][1] + "\n" + fq_reads["r1"][2] + "\n" + fq_reads["r1"][3] + "\n")
    r2file.write(fq_reads["r2"][0] + "\n" + fq_reads["r2"][1] + "\n" + fq_reads["r2"][2] + "\n" + fq_reads["r2"][3] + "\n")

# close the input files 
fq1.close()
fq2.close()
fq3.close()
fq4.close()

# close all the output FASTQ files 
for file in file_handlers_dict:
    file_handlers_dict[file].close()

# stats
total_reads = counter_matched + counter_swapped + counter_unknown
perc_matched = counter_matched / total_reads * 100
perc_swapped = counter_swapped / total_reads * 100
perc_unknown = counter_unknown / total_reads *100


with open(countsfile, "w") as fw:
    fw.write("Number (& %) of matched index-pairs: " + str(counter_matched) + " " + str(perc_matched) + "%\n")
    fw.write("Number (& %) of swapped index-pairs: " + str(counter_swapped) + " " + str(perc_swapped) + "%\n")
    fw.write("Number (& %) of unknown index-pairs: " + str(counter_unknown) + " " + str(perc_unknown) + "%\n")
    fw.write("Matched read-pairs: \n") 
    for i in counter_matched_dict:
        fw.write(i + "\t" + str(counter_matched_dict[i]) + "\n")
    fw.write("Swapped read-pairs: \n")
    for i in counter_swapped_dict:
        fw.write(i + "\t" + str(counter_swapped_dict[i]) + "\n") 
    

# percentage of reads with from each sample (A1 = #reads/ total reads)
# percentage of index hopped, matching, and swapped

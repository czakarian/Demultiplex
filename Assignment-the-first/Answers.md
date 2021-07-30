# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    ```Read 1```
    ![Read 1](https://github.com/czakarian/Demultiplex/blob/master/Assignment-the-first/R1_plot.png)
    ```Index 1```
    ![Read 2](https://github.com/czakarian/Demultiplex/blob/master/Assignment-the-first/R2_plot.png)
    ```Index 2```
    ![Read 3](https://github.com/czakarian/Demultiplex/blob/master/Assignment-the-first/R3_plot.png)
    ```Read 2```
    ![Read 4](https://github.com/czakarian/Demultiplex/blob/master/Assignment-the-first/R4_plot.png)

What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
```
The quality score cutoff for the index reads should be stricter than the quality score cutoff for biological reads
since our demultiplexing is dependent on extracting out any reads that have been index swapped or have unknown indexes.
This means it is critical to be certain that the reported index sequence is actually the expected index sequence.
A quality score of 30 means there is only a 0.01% chance that the base call is incorrect, which should be sufficient 
for index reads. For biological reads a lower quality score should suffice such as 20, which means there is a 1%
chance that the base call is incorrect.
```

How many indexes have undetermined (N) base calls?
```
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
        3976613
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
        3328051
```

    
## Part 2
1. Define the problem
```
    We want to develop a strategy to perform demultiplexing of RNA seq results 
    that were returned in the form of FASTQ files. Demultiplexing involves
    the removal of read pairs that have been index swapped or have
    unknown indexes and the organization of all the read pairs with matched indices
    into specific buckets corresponding to their index. Given 4 total FASTQ files,
    a pair of read files and a pair of index files, we want to separate
    each read pair according to which of the following categories it falls into:
    matched index pair, swapped index pair, unknown index pair. We are also
    given a text file containing the list of indexes that were included in
    the dual matched library preparation which will help you categorizes
    unknown indexes that arise.
```
2. Describe output
```
    A pair of FASTQ files for each matching index-pair, one for R1 and one for R2 [24 x 2 = 48 FASTQ files]
    A pair of FASTQ files for swapped index-pairs (1 x 2 = 2 FASTQ files)
    A pair of FASTQ files for reads with unknown or low quality indices (1 x 2 = 2 FASTQ files)
    The counts for the # read pairs with matching indices (per index-pair), read pairs with index-hopping, and read pairs with unknown indices
```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
```
    To get a list of indexes, open and parse through the indexes text file
        Store the index and the index label in a dictionary [will be referenced later when naming the output FASTQ files)
        Do we want to store the other info in each row?
    Open up all the output FASTQ files (48 matched files + 2 swapped + 2 unknown) to write to [don't use 'with open' just 'open']
        there should be R1, R2 pair for each index, for swapped, and for unknown

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
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
```
    def revcomp(seq:str) -> str
    """This function returns the reverse complement of a DNA sequence
        returns reverse complement sequence as a string"""
        # Define a dictionary and fill it with nucleotides and their complements 
            each key is a nucleotide and value is its complement nucleotide
        # Take input str and produce a new str where each nucleotide is replaced by its complement
        # Produce a new str with the complemented string reversed 
        return variable storing the reverse complemented sequence
    Example:
        Input: ATTGGC
        Expected output: GCCAAT

    convert_phred(score:str) -> int
    """This function takes a phred score and converts it to a quality score"""
        returns quality score as an integer
    # Using ord function convert ASCII value of phred score to an integer and add 33 
    Example:
        Input: I
        Expected output: 40
```

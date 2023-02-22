# Demultiplexer

Given four input fastq files (2 with biological reads, 2 with index reads) and a list of  known indexes, this program will demultiplex reads by index-pair, outputting one R1 fastq file and one R2 fastq file per matching index-pair, another two fastq files for non-matching index-pairs (index-hopping), and two additional fastq files when one or both index reads are unknown or low quality.
 
The sequence of each index-pair will be added to the header of BOTH reads in all fastq files for all categories (e.g. “AAAAAAAA-CCCCCCCC” will be appended to headers of every read pair that had an index1 of AAAAAAAA and an index2 of CCCCCCCC.

Final output stats files will report the number of read-pairs with properly matched indexes (per index-pair), the number of read pairs with index-hopping observed, and the number of read-pairs with unknown index(es). 

### Input
1. 4 fastq files (one read pair, one index pair)
2. A text file with a list of known index sequences

- argparse options:
    - ```-f```, ```--files```: required arg, Paths to input fastq files (one read pair, one index pair)
    - ```-i```, ```--indexes```: required arg, Path to file containing known index sequences + sample information
    - ```-d```, ```--direct```: required arg, Path to output directory
    - ```-s```, ```--stats```: required arg, Name for output stats files
    
### Output
1. A pair of fastq files per known index pair, a pair for index-hopped read-pairs, and a pair for reads with unknown or low quality index-pairs
2. Summary stats files
    - % and # of read-pairs with matched indexes, read-pairs with index-hopping, and read-pairs with unknown indexes
    - % and # of read-pairs with matched indexes reported per index-pair

#Bioinfo module 

DNA_bases = ["A","T","G","C","N","a","t","g","c","n"]
RNA_bases = ["A","U","G","C","N","a","u","g","c","n"]

def validate_base_seq(seq: str,RNAflag: bool =False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score given a string of quality scores"""
    qsum = 0
    for letter in phred_score:
        qsum += convert_phred(letter)
    return qsum/len(phred_score)

def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters"
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_N50(ls: list) -> int:
    """This function takes as input a list of contig lengths and returns the N50 value of the assembly"""
    midpoint = 0.5 * sum(ls) # store the midpoint of the assembly
    ls_sorted = sorted(ls,reverse=True) # sort lengths in order from largest to smallest
    s = 0 # will store the sum of lengths until s >= midpoint 
    for i in ls_sorted:  
        s += i 
        if s >= midpoint:
            return i
    return 0

def oneLineFasta(input: str, output: str) -> None:
    """This function takes a FASTA file with multiple sequence lines per read and
    returns a FASTA file with one sequence line per read."""
    with open(input, "r") as fr, open(output, "w") as fw:
        first_header = True
        for line in fr:
            # if we are at a header line, output the header to file
            if line[0:1] == ">": 
                # make sure to add a new line char before each header (expect the first) since we strip the new line from end of the prev seq line
                if first_header != True:
                    line = "\n" + line
                first_header = False         
                fw.write(line)
            # if we are at a seq line, write to the file but strip new line char
            else:
                fw.write(line.strip())


if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("validate_base_seq passed DNA and RNA tests")

    assert convert_phred("#") == 2, "convert_phred returns incorrect quality score value"
    assert convert_phred("5") == 20, "convert_phred returns incorrect quality score value"
    assert convert_phred("I") == 40, "convert_phred returns incorrect quality score value"
    print("convert_phred correctly converted phred scores")

    assert qual_score("FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@") == 37.62105263157895, "qual_score does not return correct average quality score value"
    print("qual_score correctly calculated  average quality score ")

    assert gc_content("ATGCGCGCTTAATTAA") == 0.375 , "gc_content does not return correct value"
    print("gc_content correctly calculated GC content")

    assert calc_N50([4,2,3,4,5,1,1]) == 4, "calc_N50 does not return correct N50 value"
    print("N50 calculated correctly")

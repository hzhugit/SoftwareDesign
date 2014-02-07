# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Huanzhen David Zhu
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    # Create a count for the DNA for easy segregation.
    dna_count = range(0, len(dna),3)

    # Initially declare an empty string amino that stores the translated amino acid.
    amino = ''
    
    # Removes the last few nucleotides if there's not enough for one codon.
    if (len(dna) - dna_count[-1]) % 3 != 0:
        del dna_count[-1]

    # Packages and translates aa to codons.
    for nucleotide in dna_count:

        # Initial split.
        cod = dna[nucleotide:nucleotide+3]

        # Translation part 1.
        for i in range(len(codons)):
            if cod in codons[i]:
                amino = amino + aa[i]
                break

    # Return translation.
    return amino

def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """
    inp = 'GACGTGGCACTGGCGT'
    print "input: " + inp
    print "expected output: DVALA"
    print "actual output: " + coding_strand_to_AA(inp)
# coding_strand_to_AA_unit_tests()

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    # Initial reverse.
    rev_dna_1 = dna[::-1]

    # Final translation storage.
    rev_dna_2 = ''
    
    # Translation.
    for nucleotide in rev_dna_1:
        if nucleotide == 'A':
            rev_dna_2 = rev_dna_2 + 'T'
        elif nucleotide == 'T':
            rev_dna_2 = rev_dna_2 + 'A'
        elif nucleotide == 'C':
            rev_dna_2 = rev_dna_2 + 'G'
        elif nucleotide == 'G':
            rev_dna_2 = rev_dna_2 + 'C'
    
    return rev_dna_2

def get_reverse_complement_unit_tests():
    """ Unit tests for the get_complement function """
    print "input: ATGCCCGCTTT"
    print "expected output: AAAGCGGGCAT"
    print "acutal output: " + get_reverse_complement("ATGCCCGCTTT")
    print "input: CCGCGTTCA"
    print "expected output: TGAACGCGG"
    print "actual output: " + get_reverse_complement("CCGCGTTCA") 
# get_reverse_complement_unit_tests()

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    # Define stop codons.
    stops = ['TAG', 'TAA', 'TGA']
    
    # Definte an empty array for storing stop codon positions.
    lowest = []
    
    # Find and test the earliest stop codon.
    for stop in stops:
        if stop in dna:
            if dna.index(stop) % 3 == 0:
                lowest.append(dna.index(stop))
    
    return dna[0:min(lowest)]

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
    print "input: ATGTGAA"
    print "expected output: ATG"
    print "actual output: " + rest_of_ORF("ATGTGAA")
    print "input: ATGAGATAGG"
    print "expected output: ATGAGA"
    print "acutal output: " + rest_of_ORF("ATGAGATAGG")
# rest_of_ORF_unit_tests()
        
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE        
     
def find_all_ORFs_oneframe_unit_tests():
    """ Unit tests for the find_all_ORFs_oneframe function """

    # YOUR IMPLEMENTATION HERE

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
        
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """

    # YOUR IMPLEMENTATION HERE

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""

    # YOUR IMPLEMENTATION HERE

def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """

    # YOUR IMPLEMENTATION HERE

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    # YOUR IMPLEMENTATION HERE

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """

    # YOUR IMPLEMENTATION HERE
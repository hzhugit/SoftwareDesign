# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Huanzhen David Zhu
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons
import random

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
    dna_count = range(0, len(dna), 3)

    # Initially declare an empty string amino that stores the translated amino acid.
    amino = ''
    
    # Removes the last few nucleotides if there's not enough for one codon.
    if (len(dna) - dna_count[-1]) % 3 != 0:
        del dna_count[-1]

    # Packages and translates aa to codons.
    for nucleotide in dna_count:

        # Initial split.
        triple = dna[nucleotide:nucleotide+3]

        # Translation part 1.
        for i in range(len(codons)):
            if triple in codons[i]:
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
    # Filter blocks.
    for nucleotide in range(0,len(dna),3):
        if dna[nucleotide:nucleotide+3] == 'TAG' or dna[nucleotide:nucleotide+3] == 'TAA' or dna[nucleotide:nucleotide+3] == 'TGA':
            return dna[0:nucleotide]
            break
    
    # Return whole sequence if not complete.
    return dna

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
    print "input: ATGTGAA"
    print "expected output: ATG"
    print "actual output: " + rest_of_ORF("ATGTGAA")
    print "input: ATGAGATACGATGAGATAGG"
    print "expected output: ATGAGATACGATGAGATAGG"
    print "acutal output: " + rest_of_ORF("ATGAGATACGATGAGATAGG")
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
    # Declare for multiple returns.
    return_array = []
    
    # For loop that finds all start codons and returns their end ORFs.
    for nucleotide in range(0,len(dna),3):
        if dna[nucleotide:nucleotide+3] == 'ATG':
            return_array.append(rest_of_ORF(dna[nucleotide:]))

    return return_array
            

def find_all_ORFs_oneframe_unit_tests():
    """ Unit tests for the find_all_ORFs_oneframe function """
    return find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCCATTTAA")
# print find_all_ORFs_oneframe_unit_tests()

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    # Declare empty return array.
    return_array = []

    # Loop through the three possible states.
    for shift in range(3):
        return_array.extend(find_all_ORFs_oneframe(dna[shift:]))

    return return_array

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
    return find_all_ORFs("ATGCATGAATGTAG")
# print find_all_ORFs_unit_tests()

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """
    return find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
# print find_all_ORFs_both_strands_unit_tests()

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""
    if len(find_all_ORFs_both_strands(dna)) == 0:
        return ''
    else:
        return max(find_all_ORFs_both_strands(dna), key=len)

def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """
    # return longest_ORF("ATGCGAATGTAGCATCAAA")
    return longest_ORF("12321312ATGATG")
# print longest_ORF_unit_tests()

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # Provide some storage space for the maxes.
    max_ORF = []

    # Shuffle and save all the longest ORFs.
    for _ in range(num_trials):
        dna = list(dna)
        random.shuffle(dna)
        new_dna = ''.join(dna)
        max_ORF.append(len(longest_ORF(new_dna)))

    return max(max_ORF)

# print longest_ORF_noncoding('ATGCATGAATGTAGATAGATGTGCCCATTTAA',1)

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """
    ORFs = find_all_ORFs_both_strands(dna)
    AAs = []

    for ORF in ORFs:
        if len(ORF) >= threshold:
            AA = coding_strand_to_AA(ORF)
            AAs.append(AA)
    
    return AAs
# print gene_finder('ATGCATGAATGTAGATAGATGTGCCCATTTAA', 30)

if __name__ == "__main__":
    from load import load_seq
    dna = load_seq("./data/X73525.fa")

    # Uncomment this find our threshold
    # longest = longest_ORF_noncoding(dna, 1500)
    # print longest
    
    longest = 891
    print gene_finder(dna, longest)
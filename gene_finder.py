# -*- coding: utf-8 -*-
"""
Week on of Gene_finder project. Exploeres some of the functions and doctesting.
Also looks at the start of the gene biology needed.

@author: Gretchen Rice

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'This is not possible'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    This one should be fine because it check strings of different lengths and
    checks different strings
    """
    revDna = list(dna)  # creates a list of all the terms in dna
    revDna.reverse()  # reverses the order of the list
    revComp = revDna  # makes sure new list is the right length
    for i in range(0, len(dna)):
        revComp[i] = get_complement(revDna[i])  # assigns values to new list

    return "".join(revComp)


def find_start(dna):
    """
        Takes a strand of DNA and looks for 'ATG' at the start of the string

        dna: string of DNA
        return: index of the start codon or -1 if none

        >>> find_start('ATGGGGGGGG')
        0
        >>> find_start('GGGGGGATGGGGGGGG')
        6
        >>> find_start('GGGAAAT')
        -1
    """
    for i in range(0, len(dna), 3):  # checks at every 3 for ATG
        c = dna[i:i+3]
        if c == 'ATG':
            return i  # if find "ATG" return the index

    return -1  # if no ATG return -1


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'

    This one should also work because it checks different length strings and
    tests ones with string not divisible by 3
    """
    i = 3  # know it starts with ATG so just kskip to index after that

    while i < len(dna):
        if dna[i] == 'T':  # check for first of end codon
            # check for rest of end codon and return dna strand
            if (dna[i+1] == 'A') and (dna[i+2] == 'G' or dna[i+2] == 'A'):
                return dna[:i]
            elif dna[i+1] == 'G' and dna[i+2] == 'A':
                return dna[:i]

        i += 3
    return dna  # return dna if no end codon


"""
My first try at find_all_Orfs which worked alone
Did not work with later code because my index was out of range
Lead me to make a new function to find start codon seperately

def find_all_ORFs_oneframe_try(dna):
     Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    orfs = []

    i = 0

    while i < len(dna):
        if (dna[i] == 'A') and (dna[i+1] == 'T') and (dna[i+2]):
            strand = rest_of_ORF(dna[i:])
            orfs.append(strand)
            i = i + len(strand)
        else:
            i += 3

    return orfs
"""


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    This one should also work to show enough because it has many ATGs, not all
    starting on something divisible by 3
    """
    orfs = []

    while True:
        startIndex = find_start(dna)  # finds ATG, start codon
        if startIndex == -1:
            break  # if no start codon, end
        else:
            newDna = dna[startIndex:]  # if start codon, create dna at index
            strand = rest_of_ORF(newDna)  # find the dna strand
            orfs.append(strand)  # add strand to new list
            newStart = len(strand)  # find new place to look for start codons
            dna = dna[newStart:]  # modify dna starting at new start pos

    return orfs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']

    Added second so tried one with length not divisble by 3
    """
    allOrfs = []
    allOrfs = allOrfs + find_all_ORFs_oneframe(dna)
    allOrfs = allOrfs + find_all_ORFs_oneframe(dna[1:])
    allOrfs = allOrfs + find_all_ORFs_oneframe(dna[2:])

    allOrfsNoDupes = []
    for orf in allOrfs:
        if orf not in allOrfsNoDupes:
            allOrfsNoDupes.append(orf)
    return allOrfsNoDupes


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    This seems good, has ATGs at different places and is length divisble by 3
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()

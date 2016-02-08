# -*- coding: utf-8 -*-
"""
Part One of Gene Finder Code

@author: Rachel Yang

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

    Needed to add two more doctests to ensure that all four nucleotides returned 
    the correct complement nucleotide

    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide in ('A'):
        return 'T'
    elif nucleotide in ('T'):
        return 'A'
    elif nucleotide in ('C'):
        return 'G'
    elif nucleotide in ('G'):
        return 'C'
    else:
        return None

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    These two tests are sufficient for this doctest because it parses through
    all four nucleotides and is evident whether or not the function works with 
    both reversing the string and finding the correct complement nucleotide.
    Adding any more tests would be redundant.
    """
    index = len(dna) - 1
    reverse_complement = ''
    while index >= 0:
        nucl_reverse_compl = dna[index]
        if nucl_reverse_compl in ('A'):
            new_nucl_RC = 'T'
        elif nucl_reverse_compl in ('T'):
            new_nucl_RC = 'A'
        elif nucl_reverse_compl in ('C'):
            new_nucl_RC = 'G'
        elif nucl_reverse_compl in ('G'):
            new_nucl_RC = 'C'
        reverse_complement = reverse_complement + new_nucl_RC
        index = index - 1
    return reverse_complement

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

    Needed to add a doctest for the third stop codon, 'TAA'
    >>> rest_of_ORF("ATGTAA")
    'ATG'

    Needed to add a doctest to test the case for "if there is no in frame stop codon, 
    returns the whole string"
    >>> rest_of_ORF("AGATGT")
    'AGATGT'
    """

    for i in range(0, len(dna),3):
        codon = dna[i:i+3]
        if codon in ['TGA','TAG','TAA']:
            return dna[:i]
    return dna

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
    """
    list_ORFs_one = []
    i = 3           
    while i < len(dna):
        if dna[i-3:i] == 'ATG':
            list_ORFs_one.append(rest_of_ORF(dna[i-3:]))
            i += len(rest_of_ORF(dna[i-3:]))
        else:
            i += 3
    return list_ORFs_one

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
    """
    one = find_all_ORFs_oneframe(dna[:])
    two = find_all_ORFs_oneframe(dna[1:])
    three = find_all_ORFs_oneframe(dna[2:])
    return one + two + three

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    strand_one = find_all_ORFs(dna)
    strand_two = find_all_ORFs(get_reverse_complement(dna))
    return strand_one + strand_two

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    step_one = find_all_ORFs_both_strands(dna)
    max_length = max(len(s) for s in step_one)
    step_two = [s for s in step_one if len(s) == max_length]
    return ''.join(step_two)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_ORF_list = []
    for i in range(num_trials):
        shuffled_dna = shuffle_string(dna)
        longest_shuffled = find_all_ORFs_both_strands(shuffled_dna)
        longest_specific = len(longest_shuffled)
        longest_ORF_list.append(longest_specific)
    return max(longest_ORF_list, key = lambda x: int(x))

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
    new_AA_list = []
    #for i in range(0,len(dna),3):
        #codon = dna[i:i+3]
        #amino_acid = aa_table[codon]
        #new_AA_list.append(amino_acid)
    #return ''.join(new_AA_list)

    i = 3
    while i <= len(dna):
        codon = dna[i-3:i]
        amino_acid = aa_table[codon]
        new_AA_list.append(amino_acid)
        i += 3
    return ''.join(new_AA_list)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        
        Can not test with a doctest because the shuffling will give a different
        answer everytime--not consistent operation, which means that there
        can not be an expected answer.

        Instead, test through with logic, with incremental programming, with checking with
        a friend/rubber-ducking it with a friend.
    """
    threshold = longest_ORF_noncoding(dna,1500)
    gene_finder_list = []
    for i in find_all_ORFs_both_strands(dna):
        if len(i) > threshold:
            gene_finder_list.append(coding_strand_to_AA(i))
    return gene_finder_list

if __name__ == "__main__":
    import doctest
    doctest.testmod()

from load import load_seq
dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))

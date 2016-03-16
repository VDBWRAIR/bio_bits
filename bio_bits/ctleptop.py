#!/usr/bin/env python
# encoding: utf-8
"""
Created by Dereje Jima on May 21, 2015
"""
from __future__ import division
from __future__ import print_function
from Bio.Seq import *
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
#from itertools import groupby
from Bio.Data import CodonTable
from Bio.Data.IUPACData import ambiguous_dna_values
#import yaml
import argparse
from bio_bits import degen
from functools import partial
from tabulate import tabulate
from bio_bits.compat import zip
import re
import sys

__docformat__ = "restructuredtext en"

AMBICODON = {"R": ["A", "G"], "Y": ["C", "T"],
             "W": ["A", "T"], "S": ["G", "C"],
             "K": ["T", "G"],
             "M": ["C", "A"], "D": ["A", "T", "G"],
             "V": ["A", "C", "G"], "H": ["A", "C", "T"],
             "B": ["C", "G", "T"], "N": ["A", "C", "T", "G"]}

def getNearbyChars(nt):
    """(str)->(list)
    >>>getNearbyChars("R")
    ['A', 'G']
    >>>getNearbyChars("Y")
    ['C', 'T']
    >>>getNearbyChars("A")
    ['A']
    """
    return AMBICODON.get(nt) or nt

def nearbyPermutations(letters, index=0):
    """(str)->(set)
    >>>nearbyPermutations("AAR")
    set(['AAG', 'AAA'])
    >>>nearbyPermutations("ARR")
    set(['AGG', 'AAG', 'AAA', 'AGA'])
    nearbyPermutations("AAA")
    set(['AAA'])
    """
    if (index >= len(letters)):
        return set([''])
    subWords = nearbyPermutations(letters, index + 1)
    nearbyLetters = getNearbyChars(letters[index])
    return permutations(subWords, nearbyLetters)

def permutations(subWords, nearbyLetters):
    """(set, list) -> (set)
    >>>permutations(set(['CA']), ['A', 'T'])
    set(['ACA', 'TCA'])
    """
    permutations = set()
    for subWord in subWords:
        for letter in nearbyLetters:
            permutations.add(letter + subWord)
    return permutations

def getaalist(codonlist):
    """(list) -> (list)
    Return aa list from a a given nt codon list.
    >>>getaalist(['AAA','ACT'])
    ['K', 'T']
    """
    aalist = []
    for codon in codonlist:
        aa = Seq(codon, IUPAC.unambiguous_dna)
        aa = str(translate(aa))
        aalist.append(aa)
    return aalist

def list_overlap(list1, list2):
    """(str, list) -> bool
    Return True  if the two list hava element that overlaps.

    >>>list_overlap('RAC',['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
    True
    >>>list_overlap('ACT',['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
    False

    """
    for i in list1:
        if i in list2:
            return True
    return False

def access_mixed_aa(file_name):
    """(str) ->(list,list,list,list).
    Return a list of amino acide code for ambiguous dna codon, position of
    ambiguous nt codon, aa name,seq id from fasta header  by reading multifasta
    nucleotide fasta file
    """
    from Bio import SeqIO
    aa = []
    nucleotide_idx = []
    nucl_codon = []
    seqids = []
    for seq_record in SeqIO.parse(file_name, 'fasta'):
        seq_id = seq_record.id
        seqline = str(seq_record.seq).upper()
        seqline = seqline.replace("-", "N")
        n = 3
        codon_list = dict( (i + n , seqline[i:i + n]) for i in range(0, len(seqline), n))
        ambi_nucl = AMBICODON.keys()
        for key, codon in sorted(codon_list.items()):
            if list_overlap(codon, ambi_nucl):
                d, e, f = codon
                m = [d, e, f]
                items = [i for i in m if i in ambi_nucl]
                indexm = m.index(items[0])
                for idx, val in enumerate(items):
                    codonlist = list(nearbyPermutations(codon))
                    val = getaalist(codonlist)
                    # remove if aa codon is the same eg. ['D', 'D']
                    val = set(val)
                    val = "/".join(sorted(val))   # yeild 'I/L'

                    key = key - 2 + indexm
                    if '/' in val:
                        nucleotide_idx.append(key)
                        nucl_codon.append(codon)
                        seqids.append(seq_id)
#                    if "/" in val and indexm == 2:
#                        key = key
#                        nucleotide_idx.append(key)
#                        nucl_codon.append(codon)
#                        seqids.append(seq_id)
#                    elif "/" in val and indexm == 1:
#                        key = key - 1
#                        nucleotide_idx.append(key)
#                        nucl_codon.append(codon)
#                        seqids.append(seq_id)
#                    elif "/" in val and indexm == 0:
#                        key = key - 2
#                        nucleotide_idx.append(key)
#                        nucl_codon.append(codon)
#                        seqids.append(seq_id)
#                    else:
#                        pass
                    aa.append(val)

            else:
                # print "codon3 ..." ,codon
                aa1 = Seq(codon, IUPAC.unambiguous_dna)
                aa1 = aa1.translate()
                aa1 = str(aa1)
                aa.append(aa1)
    #print aa, nucleotide_idx, nucl_codon, seqids
    return aa, nucleotide_idx, nucl_codon, seqids


def create_args():
    """
    Return command line arguments

    """
    parser = argparse.ArgumentParser(
        description='Convert inframe nucleotide \
             fasta file to protein and report mixed \
             (ambiguous codon) with its location in \
             the sequence',
        epilog = '%(prog)s -i tests/Den4_MAAPS_TestData16.fasta -o out_file.txt'
    )
    g = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument("-i", type=str, help="Nucleotide fasta file", required=True)
    parser.add_argument("-o", type=str,  help="output file name", required=True)
    g.add_argument("--gb-file", type=str,  help="genbank file name")
    g.add_argument("--gb-id", type=str,  help="genabnk accession id")
    g.add_argument("--tab-file", type=str,  help="gene tab/csv file")
    parser.add_argument('--cds', type=str, help="CDS start stop[start,stop]")
    return parser.parse_args()

def mod_entry(entry, cds=None):
    '''
    Find Gap positions and non-coding region positions
    :param entry: iterable of (seqid,nucindex,aaindex,nuclcodon,aacodon,genename)
    :cds: Gene of CDS info
    :return: entry modified to reflect gap or non-coding
    '''
    new_entry = list(entry)
    nuc_pos = entry[1]
    nt = entry[3]
    if cds:
        if cds.start >= nuc_pos or cds.end <= nuc_pos:
            new_entry[4] = 'NON-CODING'
        elif 'N' in nt:
            new_entry[4] = 'GAPFOUND'
    elif 'N' in nt:
        new_entry[4] = 'GAPFOUND'

    return tuple(new_entry)

def main():
    args = create_args()
    file_name = args.i
    outfile = args.o
        
    with open(outfile, 'w+') as outf:
        aa, nuc_idx, nucl_codon, seqids = access_mixed_aa(file_name)

        # Get Gene info

        # Remove all non-mixed positions
        amb_aa_codon = filter(lambda x: '/' in x, aa)
        # get amino acid index list
        amb_aa_indx = map(lambda x: x//3 + 1, nuc_idx)

        if args.cds:
            reference_genes, cds = degen.get_genes(args.gb_id, args.gb_file, args.tab_file)
            overlapped_genes = degen.get_degen_list_overlap(reference_genes, nuc_idx)
            mixed_positions = zip(seqids, nuc_idx, amb_aa_indx, nucl_codon, amb_aa_codon, overlapped_genes)
            cds_start, cds_end = map(int, args.cds.split(','))
            cds = degen.Gene('CDS', cds_start, cds_end)
            mixed_positions= map(lambda x: mod_entry(x, cds), mixed_positions)
            headers=[
                    'seq id', 'nt Position', 'aa position',
                    'nt composition', 'aa composition', 'gene name'
                ]

        else:
            mixed_positions = zip(seqids, nuc_idx, amb_aa_indx, nucl_codon, amb_aa_codon)
            mixed_positions = map(mod_entry, mixed_positions)
            headers=[
                    'seq id', 'nt Position', 'aa position',
                    'nt composition', 'aa composition'
                ]

        # mark gaps and non-coding positions
        outf.write(
            tabulate(
                mixed_positions,
                headers=headers
            ) + "\n"
        )

if __name__ == '__main__':
    main()

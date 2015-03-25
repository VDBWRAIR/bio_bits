# coding: utf-8

from os.path import *
import tempfile
import subprocess
import multiprocessing
import pickle
from glob import glob
from contextlib import nested
import os
import shutil

from Bio import SeqIO

def get_seqmapping(inputfasta):
    '''
    Return dictionary of {'Seq#_': 'origname'}
    '''
    # Holds our mapping that we will return
    mapping = {}

    # To help verify no duplicate ids in inputfasta
    seenids = set()

    # Generate our enumeration with new Seq#_ names
    for newname, seqrecord in get_fasta_seq_enumerate(inputfasta):
        # If already seen this input seqrec.id
        if seqrecord.id in seenids:
            raise ValueError('{0} has identical sequence names in it({1})'.format(inputfasta,newname))
        # Record the mapping
        mapping[newname] = seqrecord.id
        # Record that we have seen this input id
        seenids.add(seqrecord.id)
    
    return mapping

def get_fasta_seq_enumerate(inputfasta, start=0):
    '''
    Create a mapping dictionary that maps each seq.id of inputfasta
    to Seq#_

    Essentially enumerate(SeqIO.parse(inputfasta,'fasta')) but instead of
    returning 0...N, return Seq0_...SeqN_

    TODO:
        - Add name format to arguments to make it configurable of what can be
            returned for ids
        - Add support for different file types like fastq
    '''
    for i, seqrec in enumerate(SeqIO.parse(inputfasta, 'fasta'), start=start):
        yield ('Seq{0}_'.format(i), seqrec)

def make_renamed_phylip(inputfasta):
    '''
    Create a phylip formatted file from an input fasta file by
    sorting all sequence identifiers in the fasta file and then
    renaming each of them sequentially with Seq#_ where # is a digit
    
    Returns the mapping of {originalname: Seq#_ name}, phylipfile
    '''
    # Init var here in global scope
    tmpfasta = None
    try:
        # Get input name and extension for later
        inputpath, extension = splitext(inputfasta)
        # Create temporary fasta file to store renamed sequences in
        fd, tmpfasta = tempfile.mkstemp(suffix='.fasta')
        # Mapping of original to -> Seq#_
        mapping = {}
        # Write renamed fasta sequences to temp file
        with open(tmpfasta,'w') as fh:
            # Write in same order of original name
            for mappedname, seq in get_fasta_seq_enumerate(inputfasta):
                # Record the mapping
                if seq.id in mapping:
                    raise ValueError
                mapping[seq.id] = mappedname
                # Rename sequence id to mapping name(>original -> >Seq#_)
                seq.id = mappedname
                # Write sequence to tempfile
                fh.write(seq.format('fasta'))
        # Convert temp fasta to sequential phylip
        phylip_file = '{0}.phy'.format(inputpath)
        SeqIO.convert(tmpfasta, 'fasta', phylip_file, 'phylip')
    except ValueError as e:
        raise ValueError('{0} has identical sequence names in it({1})'.format(inputfasta,e.message))
    finally:
        # Remove tempfile
        if tmpfasta is not None and isfile(tmpfasta):
            os.unlink(tmpfasta)
    return mapping, phylip_file

def rename_sequences(filepath, mapping):
    '''
    Rename file's contents such that all occurance of mapping's
    values are replaced with mapping's keys

    mapping is {'origname':'Seq#_'}
    '''
    try:
        # Get a tempfile to put contents of renamed in
        fd, tmpfile = tempfile.mkstemp()
        # Open both files
        with nested(open(filepath),open(tmpfile,'w')) as (fhin,fhout):
            # Iterate all lines in the inputfile
            for line in fhin:
                # Iterate all mappings
                for origname, renamedname in mapping.iteritems():
                    # Replace all occurances of Seq#_ with
                    # the original name and write to the tempfile
                    line = line.replace(renamedname, origname)
                fhout.write(line)
        # replace original file with the renamed one
        shutil.copy(tmpfile, filepath)
    finally:
        # Remove tempfile
        os.unlink(tmpfile)

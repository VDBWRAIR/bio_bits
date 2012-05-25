#!/usr/bin/env python2.4

##########################################################################
##                             Sequence Concate				##
##	Author: Tyghe Vallard						##
##	Date: 5/25/2012							##
##	Version: 1.0							##
##	Description:							##
##	 This program will read the fa directory for all .fa files	##
##	 and the first file it uses it will read through all of the	##
##	 identifiers and find those identifiers in each of the other	##
##	 .fa files and then merge the cooresponding sequences from 	##
##	 those files by concating them onto the first one.		##
##									##
##	 EXAMPLE:							##
##	  file1.fa contents						##
##	   >identifier1							##
##	   AAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCCCCCAAAAAAAAAAA		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   >identifier2							##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   CCCCCCCTTTTTTTAAAAAAAGGGGGGGGCCCCCCAAAAAAATTTTTTTTT		##
##	  file2.fa contents						##
##	   >identifier1							##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   CCCCCCCTTTTTTTAAAAAAAGGGGGGGGCCCCCCAAAAAAATTTTTTTTT		##
##	   >identifier2							##
##	   AAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCCCCCAAAAAAAAAAA		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	  file3.fa contents						##
##	   >identifier1							##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   >identifier2							##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##									##
##	  OUTPUT(in the actual output the newlines on the sequences	##
##		 are removed:						##
##	   >identifier1							##
##	   AAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCCCCCAAAAAAAAAAA		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   CCCCCCCTTTTTTTAAAAAAAGGGGGGGGCCCCCCAAAAAAATTTTTTTTT		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   >identifier2							##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##	   CCCCCCCTTTTTTTAAAAAAAGGGGGGGGCCCCCCAAAAAAATTTTTTTTT		##
##	   AAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCCCCCAAAAAAAAAAA		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   CCCCCCCCCAAATTTTTTTTGGGGGGGGGCCCCCCCAAAAAAAAATTTTTT		##
##	   ACCCCCCCCCTTTTTTTTTAAAAAAAAAAACCCCCCCCGGGGGGGGGAAAA		##
##########################################################################
import os
from optparse import OptionParser,OptionGroup
import datetime
import sys

#==================== Variables =========================
# The fast files that are to be read
fasta_files = []

# The options for the script
options = None

# Date format
date_fmt = "%Y-%m-%d %H:%M:%S"

# Our argument parser
parser = OptionParser()

# Output types
NORMAL = 0
WARN = 1
INFO = 2
DEBUG = 3
debug_types = ( 'Normal', 'Warn', 'Info', 'Debug' )
#---------------------------------------------------------

#==================== Functions =========================
def set_opts( ):
    global parser
    global options

    parser.add_option( "-d", "--fasta_dir", dest="fasta_dir", help="The fasta directory of the sequences to merge" )
    parser.add_option( "--debug", action="store", dest="debug", default=0, help="Set debug level. Normal(default) = 0, Warn = 1, Info = 2, Debug = 3" )
    options,args = parser.parse_args()
    if not options.fasta_dir:
        parser.print_help()
        parser.error( "Need to specify the fasta directory" )

def log( msg, type=NORMAL ):
    global date_fmt
    global options
    global debug_types
    nowdate = datetime.datetime.now().strftime( date_fmt )
    if type <= int( options.debug ):
        print "[%s] [%s] %s" % (nowdate, debug_types[type], msg)

def read_fasta_dir( dir_path ):
    """ Gathers the file names of the fasta files in a directory """
    valid_ext = ('.fna', '.fas')
    fasta_files = []
    if not os.path.isdir( dir_path ):
        print( "Invalid directory given(Not found) %s" % dir_path )
        sys.exit( -1 )
    for file in os.listdir( dir_path ):
        name,ext = os.path.splitext( file )
        if ext in valid_ext:
            fasta_files.append( file )
    fasta_files.sort()
    log( "%s" % fasta_files, DEBUG )
    return fasta_files

def read_fasta_file( fasta_dir, fasta_file ):
    """ Reads a fasta file and outputs a dictionary where the key is the sequence name and the value is the sequence """
    genes = {}
    fasta_file_path = os.path.join( fasta_dir, fasta_file )
    lines = []
    if not os.path.isfile( fasta_file_path ):
        print( "Invalid file path given(Not found) %s" % fasta_file_path )
        sys.exit( -1 )
    fh = open( fasta_file_path, 'r' )
    name_line = True
    last_name = ''
    for line in fh.readlines():
        log( line, DEBUG )
        if line[0] == '>':
            log( "Setting name: %s" % last_name, DEBUG )
            last_name = line[1:].strip()
            name_line = True
            genes[last_name] = ''
        else:
            log( "Appending sequence to %s" % last_name, DEBUG )
            genes[last_name] += line.strip()
    return genes

def all_files( dir_path ):
    """
        Given a fasta directory read all the files in that directory into one big dictionary keyed by the file name and the values are
        what is returened by read_fasta_file
    """
    all_files = {}
    fasta_files = read_fasta_dir( dir_path )
    for file in fasta_files:
        all_files[file] = read_fasta_file( dir_path, file )
    return all_files

def pretty_print_files( files_dict ):
    """ Prints a fasta files dictionary in easier to read format """
    for file,genes in files_dict.iteritems():
        print( file )
        for name,seq in genes.iteritems():
            print( "-- %s => %s" % (name,seq) )

def merge_files( files_dict ):
    # Temp merged_files array
    merged_files = {}
    # Go through all the sequences in the first file and merge each of the sequences in the other files that have the same sequence names
    file_names = files_dict.keys()
    file_names.sort()
    log( "File names to merge %s" % file_names, DEBUG )
    first_file_name = file_names[0]
    first_file = files_dict[first_file_name]
    del files_dict[first_file_name]
    for name,seq in first_file.iteritems():
        log( "Setting %s from file %s to be merged" % (name,first_file_name), DEBUG )
        merged_files[name] = seq
        last_name = name
        missing = False
        # Go through each of the other files [1:] and append each sequence to it's cooresponding sequence(merge) in order
        for file_name in sorted( files_dict.keys() ):
            # The gene name
            gene = files_dict[file_name]
            # First check to see if the name exists in this file
            #  Append if it does exist
            #  Log if it does not
            if name in gene:
                log( "Merging %s from %s" % (name,file_name), DEBUG )
                merged_files[name] += seq
            else:
                log( "%s not fund in %s" % (name, file_name), WARN )
                last_name = name
                missing = True
        if missing:
            log( "Removing %s from merged list because it is missing in some files" % name, INFO )
            del merged_files[name]
    return merged_files

def merged_fasta( merged_files ):
    """ Print out the fasta format for the merged files """
    ret_string = ""
    for name,seq in merged_files.iteritems():
        ret_string = ">%s\r\n%s" % (name,seq)
    return ret_string
        
#---------------------------------------------------------
# Set the arguments of the script
set_opts()
files = all_files( options.fasta_dir )
merged_files = merge_files( files )

print merged_fasta( merged_files )

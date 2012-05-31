#!/usr/bin/env python2.4

##########################################################################
##                       Sequence Concat
##	Author: Tyghe Vallard						                        
##	Date: 5/30/2012							                            
##	Version: 1.0							                            
##	Description:							                            
##      This script will merge all sequences with identical names in a
##      file.
##
##  Example:
##      SampleFile1.fasta contents
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##      >identifier1
##      TTTTTTTTTTTTTTTTTTTTTTTTT
##      >identifier1
##      GGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCC
##      >identifier2
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##
##      OUTPUT
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAA
##      
##########################################################################
import os
from optparse import OptionParser,OptionGroup
import sys
import cStringIO

try:
    from pyWrairLib.parser import fasta
except ImportError:
    print "Please ensure you have set your python path to include the underlying directory for pyWrairLib"
    print "Probably export PYTHONPATH=$PYTHONPATH:/home/EIDRUdata/Tyghe/Repos"
    print "Exiting..."
    sys.exit( -1 )

class SequenceConcat:
    _fasta_file_path = None
    _strip_chars = ""
    _parsed_fasta = None

    def __init__( self, fasta_file, strip_chars = "-" ):
        self._fasta_file_path = fasta_file
        self._strip_chars = strip_chars
        self._read_genbank()

    def _read_genbank( self ):
        """ Reads the genbank file into a easy to work with dictionary """
        self._parsed_fasta = fasta.read_genbank_fasta( self._fasta_file_path, self._strip_chars )

    def prune_missing( self, fasta_dict, segments_expected = 8 ):
        """
            Prunes the dictionary of sequences so that only the sequences that have segments_expected
            amount of segments are included.
            
            Note:
                This is an inplace operation
            
            +++ Unit Test +++
            # Try an inhouse easy test
            >>> s = SequenceConcat( 'example_files/example1.txt' )
            >>> fasta = s.get_sequences()
            >>> pruned_fasta = s.prune_missing( fasta, 3 )
            >>> print pruned_fasta
            {'ident4': {'1': 'AAAAAAAAAAAAAAAAA', '3': 'TTTTTTTTTTTTTTTTT'}, 'ident5': {'1': 'AAAAAAAAAAAAAAAAA'}, 'ident3': {'1': 'AAAAAAAAAAAAAAAAA', '2': 'CCCCCCCCCCCCCCCCC'}}
            >>> fasta != pruned_fasta
            True
            >>> len( fasta ) == 2
            True
        """
        # We will return a dictionary that contains the sequences that have missing segments
        segments_missing = {}
        # Delete any of the sequences names from the dictionary that
        #  do not have the required amount of segments
        for seq_name in fasta_dict.keys():
            if len( fasta_dict[seq_name] ) != segments_expected:
                # Copy the sequence
                segments_missing[seq_name] = fasta_dict[seq_name]
                # Delete it
                del fasta_dict[seq_name]
        return segments_missing

    def get_sequences( self ):
        """
            Return unmodified fasta dictionary
            
            +++ Unit Tests +++
            >>> s = SequenceConcat( 'example_files/example1.txt' )
            >>> len( s.get_sequences() ) > 1
            True
        """
        return self._parsed_fasta

    def get_merged_sequences( self, prune = True, segments_required = 8 ):
        """
            Returns a merged fasta formatted file
            Set prune to false if you don't want to prune out the sequences that don't have
            the correct amount of segments
        
            +++ Unit Tests +++
            >>> s = SequenceConcat( 'example_files/example1.txt' )
            >>> print s.get_merged_sequences( True, 3 )
            >ident1
            AAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTT
            >ident2
            AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC
            <BLANKLINE>

            # Verification using sample data
            >>> s = SequenceConcat( 'example_files/WR2848N.fasta' )
            >>> fh = open( 'example_files/WR2848N_merged.fasta' )
            >>> WR2848N_manually_merged = fh.read()
            >>> fh.close()
            >>> s.get_merged_sequences( ) == WR2848N_manually_merged[:-1]
            True
        """

        # The fasta sequences
        fasta_dict = self.get_sequences()

        # Fast String writer
        cs = cStringIO.StringIO()

        # Return variable
        output = None

        # Will hold missing sequences
        missing = None

        # If the prune option is set then prune the sequences
        if prune:
            try:
                segments_required = int( segments_required )
            except ValueError:
                print( "Invalid value for required segments" )
                sys.exit( -1 )
            missing = self.prune_missing( fasta_dict, segments_required )

            if len( missing ):
                sys.stderr.write( "==================== Sequences Missing Segments  ====================\n" )

            #segments_required
            # Loop through the missing segment sequences
            for name, segments in missing.items():
                missing_segs = [str(i) for i in range( 1, segments_required + 1 ) if str(i) not in segments]
                sys.stderr.write( ">%s is missing segment[s] %s\n" % (name, ",".join( missing_segs ) ) )

        # Loop through each sequence and merge the segments
        for name, segments in fasta_dict.items():
            cs.write( fasta.merge_segments( name, segments ) )

        output = cs.getvalue()
        cs.close()
        
        return output
                
################ Script Functions #######################
def set_opts( parser ):
    """ Set script options """
    parser.add_option( "-f", "--fasta", dest="fasta_file", help="The fasta file of the sequences to merge" )
    parser.add_option( "--strip", dest="strip_chars", default="-", help="List of characters to strip from the sequences. Default is none" )
    parser.add_option( "--test", dest="test", action="store_true", help="Run the tests for this script" )
    prune_group = OptionGroup( parser, "Pruning Options" )
    prune_group.add_option( "--noprune", dest="prune", action="store_false", default=True, help="Don't prune out sequences that don't have the required amount of segments" )
    prune_group.add_option( "--segments_required", dest="segments_required", default=8, help="Required amount of segments per sequence name" )
    parser.add_option_group( prune_group )
    options,args = parser.parse_args()
    if not options.fasta_file and not options.test:
        parser.print_help()
        parser.error( "Need to specify the fasta file" )
    return options

# Run The script if this script is executed
if __name__ == '__main__':
    parser = OptionParser()
    options = set_opts( parser )
    if options.test:
        import doctest
        doctest.testmod()
    else:
        sc = SequenceConcat( options.fasta_file, options.strip_chars )
        print sc.get_merged_sequences( options.prune, options.segments_required )

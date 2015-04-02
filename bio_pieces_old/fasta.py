import os
import sys
import re
import cStringIO

class UnknownIdentifierLineException(Exception): pass

def merge_segments( sequence_name, sequence_dict, segments_expected = [1,2,3,4,5,6,7,8] ):
    """
        Merges a sequence segmented dictionary

        +++ Unit Tests +++
        >>> print merge_segments( 'ident1', { '2': 'AAAAAAAAAAAAAAAA', '3': 'TTTTTTTTTTTTTTTT', '1': 'GGGGGGGGGGGGGGG' }, range( 1, 4 ) )
        >ident1
        GGGGGGGGGGGGGGGAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT
        <BLANKLINE>
        >>> print merge_segments( 'ident1', { 'B': 'AAAAAAAAAAAAAAAA', 'C': 'TTTTTTTTTTTTTTTT', 'A': 'GGGGGGGGGGGGGGG' }, ['A','B','C'] )
        >ident1
        GGGGGGGGGGGGGGGAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT
        <BLANKLINE>
    """
    # Start a fast stringio operator
    fasta_output = cStringIO.StringIO()
    fasta_string = None

    # First write the identifier line
    fasta_output.write( ">%s\n" % sequence_name )

    # Loop through each of the segments in order
    for seg in segments_expected:
        # In case no pruning was done, we need to skip missing segments
        if str( seg ) in sequence_dict:
            fasta_output.write( sequence_dict[str(seg)] )
    
    # Newline after the sequence is finished
    fasta_output.write( "\n" )
    
    # Get the string value and close the cStringIO object
    fasta_string = fasta_output.getvalue()
    fasta_output.close()

    return fasta_string

def is_identifier_line( fasta_line ):
    """ 
        Identifier lines in fasta files begin with >

        ++ Unit Tests ++
        >>> idents = ['>goodident','badident']
        >>> for ident in idents:
        ...     is_identifier_line( ident )
        ...
        True
        False
    """
    return fasta_line[0] == '>'

def parse_genbank_identifier( identifier ):
    """
        Parse genbank identifier line 
        Format expected
        >A/Ankara/WRAIR1425T/2009 8 (NS)

        Returns a dictionary 

        ++ Unit Tests ++
        # Check to make sure that the UnknownIdentifierLineException is thrown on a bad identifier
        >>> parse_genbank_identifier( "bad_identifier_line bad" )
        Traceback (most recent call last):
          File "/usr/lib64/python2.4/doctest.py", line 1243, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.parse_genbank_identifier[0]>", line 1, in ?
            parse_genbank_identifier( "bad_identifier_line bad" )
          File "fasta.py", line 44, in parse_genbank_identifier
            raise UnknownIdentifierLineException( identifier )
        UnknownIdentifierLineException: bad_identifier_line bad

        # Check to make sure a valid identifier is parsed
        >>> print parse_genbank_identifier( ">A/Addis Ababa/WR2848N/2009 1 (PB2)" )
        {'num': '1', 'name': 'A/Addis Ababa/WR2848N/2009', 'gene': '(PB2)'}

        # Check all possibilities
        >>> print parse_genbank_identifier( ">aAzZ-.    09_/ 52 (PB2)" )
        {'num': '52', 'name': 'aAzZ-.    09_/', 'gene': '(PB2)'}

        >>> print parse_genbank_identifier( ">BRD40616N.1" )
        {'num': '', 'name': 'BRD40616N.1', 'gene': ''}
    """
    ident_pattern = ">(?P<name>[-./a-zA-Z0-9_\s]+)\s+(?P<num>\d+)\s+(?P<gene>\(\w+\))"
    cPattern = re.compile( ident_pattern )
    try:
        m = cPattern.match( identifier )
        return m.groupdict()
    except AttributeError, e:
        if " " in identifier.strip():
            raise UnknownIdentifierLineException( identifier )
        else:
            return {'num': '', 'name': identifier[1:], 'gene': ''}

def parse_gisaid_identifier( identifier ):
    """
        Parse gisaid identifier line 
        Format expected
            Example Identifier:
                identifier------------------   gisaid id-----   col date--  si sl  seg  accession#
                >A/Hubei-Wuchang/SWL174/2010 | EPI_ISL_75881 | 2010-04-11 |  |  | NP | EPI266208
            Expected Input File Format:
            identifier_line                     := <identifier> <metadata>
            carrot                              := ">"
            identifier                          := <carrot>[-./a-zA-Z0-9_\s]
            metadata                            := <divider> <meta>
            divider                             := "|"
            meta                                := [<isolate-gisaid-id> | <collection-date> | <sample-id-provider> | <sample-id-submitting-laboratory> | <segment> | <dna-accession-number>]
            isolate-gisaid-id                    := \w+
            collection-date                     := \d{4}-\d{2}-\d{2}
            sample-id-provider                  := [ \s | [-./a-zA-Z0-9_\s] ]
            sample-id-submitting-laboratory     := \d+
            segment                             := \w+
            dna-accession-number                := \w+

        Returns a dictionary keyed with meta data values

        ++ Unit Tests ++
        # Check to make sure that the UnknownIdentifierLineException is thrown on a bad identifier
        >>> parse_genbank_identifier( "bad_identifier_line bad" )
        Traceback (most recent call last):
          File "/usr/lib64/python2.4/doctest.py", line 1243, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.parse_genbank_identifier[0]>", line 1, in ?
            parse_genbank_identifier( "bad_identifier_line bad" )
          File "fasta.py", line 44, in parse_genbank_identifier
            raise UnknownIdentifierLineException( identifier )
        UnknownIdentifierLineException: bad_identifier_line bad

        # Check to make sure a valid identifier is parsed
        >>> required_keys = ['accession_number', 'collection_date', 'gisaid_id', 'name', 'sample_provider', 'sample_submitting_laboratory', 'segment']
        >>> a = parse_gisaid_identifier( ">A/Hubei-Wuchang/SWL174/2010 | EPI_ISL_75881 | 2010-04-11 |  |  | NP | EPI266208" )
        >>> sorted( a.keys() ) == required_keys
        True
        >>> print sorted( a.values() )
        ['', '', '2010-04-11', 'A/Hubei-Wuchang/SWL174/2010', 'EPI266208', 'EPI_ISL_75881', 'NP']
        >>> a =  parse_gisaid_identifier( ">A/Pakistan/154/2009 | EPI_ISL_104874 | 2009-12-02 | L-154EIC/09 ORIGINAL | 2010705413 | HA | EPI355228" )
        >>> sorted( a.keys() ) == required_keys
        True
        >>> print sorted( a.values() )
        ['2009-12-02', '2010705413', 'A/Pakistan/154/2009', 'EPI355228', 'EPI_ISL_104874', 'HA', 'L-154EIC/09 ORIGINAL']
        >>> a = parse_gisaid_identifier( ">1 | 2 | 3 | 4 | 5 | 6 | 7" )
        >>> sorted( a.keys() ) == required_keys
        True
        >>> print sorted( a.values() )
        ['1', '2', '3', '4', '5', '6', '7']
        >>> a = parse_gisaid_identifier( ">-azA /Z09_ | a_Z_0_9 | 1111-22-22 |  |  | A_Z_a_z_0_9 | a_z_A_Z_0_9" )
        >>> sorted( a.keys() ) == required_keys
        True
        >>> print sorted( a.values() )
        ['', '', '-azA /Z09_', '1111-22-22', 'A_Z_a_z_0_9', 'a_Z_0_9', 'a_z_A_Z_0_9']
        >>> a = parse_gisaid_identifier( ">-azA /Z09_ | a_Z_0_9 | 1111-22-22 | -azA /Z09_ | 012349 | A_Z_a_z_0_9 | a_z_A_Z_0_9" )
        >>> sorted( a.keys() ) == required_keys
        True
        >>> print sorted( a.values() )
        ['-azA /Z09_', '-azA /Z09_', '012349', '1111-22-22', 'A_Z_a_z_0_9', 'a_Z_0_9', 'a_z_A_Z_0_9']
        >>> a = parse_gisaid_identifier( ">-azA /Z09_ | a_Z_0_9 | 1111-22-22 |  |  | A_Z_a_z_0_9 | a_z_A_Z_0_9" )
        >>> sorted( a.keys() ) == required_keys
        True
    """
    ident_pattern = ">(?P<name>.*)\s\|\s(?P<gisaid_id>.*)\s\|\s(?P<collection_date>.*)\s\|\s(?P<sample_provider>.*)\s\|\s(?P<sample_submitting_laboratory>.*)\s\|\s(?P<segment>.*)\s\|\s(?P<accession_number>.*)"
    cPattern = re.compile( ident_pattern )
    try:
        m = cPattern.match( identifier )
        return m.groupdict()
    except AttributeError, e:
        raise UnknownIdentifierLineException( identifier )

def read_gisaid_fasta( fasta_file, strip_these="-" ):
    """ 
        Reads gisaid fasta file which contains many of the same identifiers split into many entries
            
        Return:
            A dictionary keyed with the identifier name and the value is a list of the sequences in order
            that they were read from the file

        +++ Unit Tests +++
        # Do a single gene test
        >>> path = os.path.dirname( __file__ )
        >>> a = read_gisaid_fasta( os.path.join( path, "example_files/gisaid_example1.txt" ) )
        >>> print len( a.keys() )
        2
        >>> print a.keys()
        ['ide/ a-b/asdf', 'ident1']
        >>> a = read_gisaid_fasta( os.path.join( path, "example_files/gisaid_example2.txt" ) )
        >>> print len( a.keys() )
        2
        >>> print len( a['A/Hubei-Wuchang/SWL174/2010'] )
        8
    """
    fasta = {}
    fasta_contents = None
    fh = open( fasta_file, 'r' )
    # Holds the last identifier seen
    last_ident = ""

    # Started read sequence flag
    start_read_seq = False

    for line in fh:
        # Identifier line
        if is_identifier_line( line ):
            # Every time a sequence line is encountered we should reset the start_read_seq flag
            # so we can keep track of which segment we are reading
            start_read_seq = False

            # Parse the identifier line(Throws exception if there is a malformed identifier)
            ident = parse_gisaid_identifier( line )

            # If the identifier_name has not been seen yet start a new entry
            #  and set the last_ident
            # We do this because of the multiple segments
            last_ident = ident
            if not ident['name'] in fasta:
                fasta[last_ident['name']] = {}
        else:
            if not start_read_seq:
                fasta[last_ident['name']][last_ident['segment']] = line.strip()
                start_read_seq = True
            else:
                fasta[last_ident['name']][last_ident['segment']] += line.strip()
            
    fh.close()
    return fasta

def read_genbank_fasta( fasta_file, strip_these="-" ):
    """ 
        Reads genbank fasta file which contains many of the same identifiers split into many entries
        >ident1 1 (gene1)
        AAAA
        >ident1 2 (gene2)
        GGGG
        >ident1 3 (gene3)
        CCCC
        >ident2 1 (gene5)
        TTTT
        >ident2 2 (gene6)
        CCCC
        
        Returns a dictionary keyed with the identifier name and the value is a list of the sequences in order
            that they were read from the file

        +++ Unit Tests +++
        # Do a single gene test
        >>> path = os.path.dirname( __file__ )
        >>> a = read_genbank_fasta( os.path.join( path, "example_files/genbank_example1.txt" ) )
        >>> print len( a.keys() )
        1
        >>> print len( a['A/Addis Ababa/WR2848N/2009'] )
        8
        >>> a = read_genbank_fasta( os.path.join( path, "example_files/genbank_example2.txt" ) )
        >>> print len( a.keys() )
        192
        >>> print len( a )
        192
    """
    fasta = {}
    fasta_contents = None
    fh = open( fasta_file, 'r' )
    # Holds the last identifier seen
    last_ident = ""

    # Started read sequence flag
    start_read_seq = False

    for line in fh:
        # Identifier line
        if is_identifier_line( line ):
            # Every time a sequence line is encountered we should reset the start_read_seq flag
            # so we can keep track of which segment we are reading
            start_read_seq = False

            # Parse the identifier line(Throws exception if there is a malformed identifier)
            ident = parse_genbank_identifier( line )

            # If the identifier_name has not been seen yet start a new entry
            #  and set the last_ident
            # We do this because of the multiple segments
            last_ident = ident
            if not ident['name'] in fasta:
                fasta[last_ident['name']] = {}
        else:
            if not start_read_seq:
                fasta[last_ident['name']][last_ident['num']] = line.strip()
                start_read_seq = True
            else:
                fasta[last_ident['name']][last_ident['num']] += line.strip()
            
    fh.close()
    return fasta

def read_fasta_file( fasta_file_path, strip_these="-" ):
    """ Reads a fasta file and outputs a dictionary where the key is the sequence name and the value is the sequence """
    genes = {}
    lines = []
    #if not os.path.isfile( fasta_file_path ):
    #    print( "Invalid file path given(Not found) %s" % fasta_file_path )
    #    sys.exit( -1 )
    fh = open( fasta_file_path, 'r' )
    name_line = True
    last_name = ''
    for line in fh.readlines():
        #log( "Line being read: %s" % line, DEBUG )
        if line[0] == '>':
            last_name = fix_name( line[1:] )
            #log( "Starting new sequence name: %s" % last_name, DEBUG )
            name_line = True
            genes[last_name] = ''
        else:
            #log( "Appending sequence %s to %s" % (line, last_name), DEBUG )
            genes[last_name] += strip_chars( line.strip(), strip_these )
    return genes

def strip_chars( sequence, strip_chars ):
    """ Strip all occurances of strip_chars from sequence """
    pattern = "|".join( [c for c in strip_chars] )
    compiled_pattern = re.compile( pattern )
    return compiled_pattern.sub( '', sequence )

def fix_name( seq_name ):
    """ Remove everything after a period and change - to _ """
    seq_name_fixed = re.sub( '-', '_', seq_name.strip() )
    return re.sub( '\..*', '', seq_name_fixed )

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    # Import the path just below the current script's path
    sys.path.append( '../' )
    _test()

'''
Tools for geting info about and comparing Variant Calling Format files.

Usage:
    vcf_utils view 
    vcf_utils filter [<FIELD> (--minimum <MIN> | --maximum <MAX> | --exists | --not-exists | --ambiguous | --not-ambiguous)]... <FILE>
    . filter_vcf [CBD --minimum MIN]...  <FILE> (#allow no filter)
    . filter_vcf [ALT --exists]
    . diff [<FIELD> --threshold <THRESH>] <FILE1> <FLIE2>
    . mutations  <FILE> (# this creates a sub-df which can then be diffed.)

Arguments:
    <FIELD>   CB, CBD, etc.
   
Options:
    -c, --count   Get the lenght of output (number of records) instead of records
    --threshold <THRESH>  Minimum difference for something to show up in diff output
    --threshold-percent <TPCNT>  THRESH as percent rather than raw value
    --minimum <MIN>   The minimum value of a given field for it to show up in output
    --maximum <MAX>   See --minimum
    --exists       Filter on presence (not '-')
    --not-exists   Opposite of --exists
    -o <DIR>, --outdir <DIR>      Ouptut to this directory [default: .]
    --ambiguous   Filter on the field is an ambiguous call (Y, etc.)
    --not-ambiguous  Opposite of ambiguous 
''' 


def get_function(string):
    getattr(vcf_compare, string)

def parse_args(docargs):
    pass

#TODO: Need default VCF headers.
#TODO: Does expected output include headers? this will screw up --count
#TODO: use PyVcf vcf_filter?

'''
If the type of a field is int/float, use an int-based threshold.
If the type of a field is otherwise (list/str), check equality.
'''

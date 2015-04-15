from schema import Schema, Use, And
from docopt import docopt
import vcfcat




def validate_args():
    '''
    Have to handle transforming ints to ints;
    '''
    pass

def main():
    raw_args = docopt(vcfcat.__doc__, version='Version 0')
    parsed_args = validate_args(raw_args)
    #TODO: dispatch here
    #return amos2fastq.make_fastqs_by_contigs(parsed_args['<fastqs>'], parsed_args['--amos'])


if __name__ == "__main__":
    main()

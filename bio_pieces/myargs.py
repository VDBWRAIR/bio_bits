'''
Usage: amos2fastq <fastqs>... --amos=<amos> 

Given any number of fastq files and a specified AMOS file (usually with .afg extension) describing a series of contigs, create seperate fastq files arranged by the contig each read aligned to.

'''
from schema import Schema, Use, And, SchemaError
from docopt import docopt
import os
import amos2fastq
#Do file validation immediately when script is started
def is_fastq(filename):
    pass

def is_amos(filename):
    pass 

def all_elemnts_unique(collection):
    return len(collection) == len(set(collection))

# possibly try/catch
def validate_args(raw_args): 
    s = Schema({'<fastqs>' : And([os.path.exists],all_elemnts_unique, error='fastq files must be different and exist'),
        '--amos' : Use(open, error='AMOS file should be readable') })
    return s.validate(raw_args)

def main():
    raw_args = docopt(__doc__, version='Version 0')
    parsed_args = validate_args(raw_args)
    return amos2fastq.main(parsed_args['<fastqs>'], parsed_args['--amos'])


if __name__ == "__main__": 
    main()



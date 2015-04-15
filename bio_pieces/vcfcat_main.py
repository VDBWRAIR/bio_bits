'''
Usage:
    vcfcat filter <FILE1> [ --tag=<TAG> (--ge | --le | --geq | --leq | --eq | --neq) <VALUE> ] [-c]
    vcfcat exists <FILE1> (--tag=TAG> )
    vcfcat ambiguous <FILE1>
    vcfcat diff <FILE1> <FILE2> (--tag=<TAG> [--threshold=<THRESH>])  [-c]
    vcfcat stat <FILE1>
    vcfcat statdiff <FILE1> <FILE2>

Options:
    -t=<THRESH> --threshold=<THRESH>  Numerical threshold to determine difference [Default: 0]
    -c, --count  Output the number of records.
    --no-header  Output without the VCF header info.
'''

from schema import Schema, Use, And, Optional
from docopt import docopt
import operator
import vcf_compare as vcfcat
import sys
import vcf
ops = ['--ge', '--le', '--geq' , '--leq' , '--eq' , '--neq']
#TODO: Find a better way to dispatch commandline apps
def validate_vcf_file(name):
    if name is None:
        return None
    f  = open(name)
    return list(vcf.Reader(f))
#    return And( Use(open, error='could not open vcf file,'),
#               Use(vcf.Reader, error='pyvcf could not read file'),
#               Use(list))

def run(raw):
    commands = ['filter', 'diff', 'stat', 'statdiff', 'exists', 'ambiguous']
    schema_dict=    {
            '<FILE1>' : Use(validate_vcf_file),
            Optional('<FILE2>') : Use(validate_vcf_file),
            '<VALUE>' : lambda v: (v) if v.isalpha() else int(v),
            '--count' : bool,
            '--threshold' : Use(int, error='Threshold must be integer'),
            '--tag' : str
             #tags.__contains__ #    error='Tag was not valid, should be one of {0}'.format(' '.join(tags)))
        }
    schema_dict.update( dict( (arg, bool) for arg in commands + ops))
    _args = Schema(schema_dict).validate(raw)
    cmd_str = [k for (k, arg) in _args.items() if k in commands and arg][0]
    return cmd_str, _args

def dispatch_cmd(cmd_str, args):
    cmd = getattr(vcfcat, cmd_str)
    if cmd_str == 'filter':
       op = [k for k, arg in args.items() if k in ops and arg ][0]
       opfunc = getattr(operator, op[2:])
       result = vcfcat.filter(args['<FILE1>'], args['--tag'], opfunc, args['<VALUE>'])
    elif cmd_str in ['stat', 'ambiguous']:
        result = cmd(args['<FILE1>'])
    elif cmd_str == 'diff':
        result = vcfcat.diff(args['<FILE1>'], args['<FILE2>'], args['--tag'], args['--threshold'])
    else: #'statdiff'
        result = vcfcat.statdiff(args['<FILE1>'], args['<FILE2>'])
    return result

def compute():
    raw_args = docopt(__doc__, version='Version 0')
    cmd_str, args = run(raw_args)
    result = dispatch_cmd(cmd_str, args)
    return result, args

def main():
    result, args = compute()
    stdout = sys.stdout
    if args['--count']:
        stdout.write(str(len(result)))
    else:
        for rec in result:
            vcf.Writer(stdout, args['<FILE1>']).write_record(rec)

    #TODO: dispatch here
    #TODO: handle -c
    #TODO: Make sure can accept input from stdin
    #return amos2fastq.make_fastqs_by_contigs(parsed_args['<fastqs>'], parsed_args['--amos'])


if __name__ == "__main__":
    main()

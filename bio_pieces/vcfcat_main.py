'''
Usage:
    vcfcat filter <FILE1> [ --tag=<TAG> (--ge | --le | --geq | --leq | --eq | --neq) <VALUE> ] [-c]
    vcfcat exists <FILE1> (--tag=<TAG> ) [-c]
    vcfcat ambiguous <FILE1>  [-c]
    vcfcat vcall <FILE1>  [--csv | -c ]
    vcfcat diff <FILE1> <FILE2> (--tag=<TAG> [--threshold=<THRESH>])  [-c]
    vcfcat stat <FILE1>
    vcfcat statdiff <FILE1> <FILE2>

Options:
    -t=<THRESH> --threshold=<THRESH>  Numerical threshold to determine difference [Default: 0]
    -c, --count  Output the number of records.
    --no-header  Output without the VCF header info.
    --csv         Output simplified format
'''

from schema import Schema, Use,  Optional
from docopt import docopt
import operator
from bio_pieces import vcfcat
import sys
import vcf
#TODO: figure out diff output
#TODO: Find a better way to dispatch commandline apps
ops = ['--ge', '--le', '--geq' , '--leq' , '--eq' , '--neq']
def validate_value(val):
    if val is None:
        return val
    return val if val.isalpha() else int(val)

def validate_vcf_file(name):
    if name is None:
        return None
    f  = open(name)
    return list(vcf.Reader(f))
#    return And( Use(open, error='could not open vcf file,'),
#               Use(vcf.Reader, error='pyvcf could not read file'),
#               Use(list))

def run(raw):
    commands = ['vcall', 'filter', 'diff', 'stat', 'statdiff', 'exists', 'ambiguous']
    schema_dict=    {
            '<FILE1>' : Use(validate_vcf_file),
            Optional('<FILE2>') : Use(validate_vcf_file),
            Optional('<VALUE>') : Use(validate_value),
            '--count' : bool,
            '--threshold' : Use(int, error='Threshold must be integer'),
            Optional('--tag') : lambda t: True,
            '--csv' : bool
             #tags.__contains__ #    error='Tag was not valid, should be one of {0}'.format(' '.join(tags)))
        }
    schema_dict.update( dict( (arg, bool) for arg in commands + ops))
    _args = Schema(schema_dict).validate(raw)
    cmd_str = [k for (k, arg) in _args.items() if k in commands and arg][0]
    filename = raw['<FILE1>']
    return cmd_str, _args, filename

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
    elif cmd_str == 'vcall':
        result = cmd(args['<FILE1>'])
    elif cmd_str == 'exists':
        result = cmd(args['<FILE1>'], args['--tag'])
    else: #'statdiff'
        result = vcfcat.statdiff(args['<FILE1>'], args['<FILE2>'])
    return result

def compute():
    raw_args = docopt(__doc__, version='Version 0')
    cmd_str, args, filename = run(raw_args)
    result = dispatch_cmd(cmd_str, args)
    return result, args, cmd_str, filename


def print_variant_call(rec):
    print("{0}\t{1}\t{2}\t{3}".format(
        rec.CHROM, rec.POS, rec.REF, rec.INFO['CB']
    ))

def main():
    result, args, cmd, filename = compute()
    template = vcf.Reader(open(filename))
    stdout = sys.stdout
    if args['--count']:
        print(str(len(result)))
    elif args['--csv']:
        print('Reference\tPosition\tReference Base\tCalled Base')
        for rec in result:
            print_variant_call(rec)
    elif cmd.startswith('stat'):
        print(str(result))
    else:
        writer = vcf.Writer(stdout, template)
        for rec in result:
            writer.write_record(rec)

    #TODO: Make sure can accept input from stdin

if __name__ == "__main__":
    main()

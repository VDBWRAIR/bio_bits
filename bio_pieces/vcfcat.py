from __future__ import print_function
'''
Usage:
    vcfcat filter <FILE> [ --tag  <TAG> (--ge | --le | --geq | --leq | --eq | --neq) <VALUE> ]... [-c]
    vcfcat filter <FILE> (--tag <TAG> --exists | --ambiguous)
    vcfcat diff <FILE1> <FILE2> (--tag <TAG> [--threshold <THRESH>])  [-c]
    vcfcat stat <FILE>
    vcfcat statdiff <FILE1> <FILE2>

Options:
    -t THRESHOLD --threshold=THRESHOLD  Numerical threshold to determine difference [Default: 0]
    -c, --count  Output the number of records.
    --no-header  Output without the VCF header info.
    --

    -o DIR --output=DIR  The output directory [default: vcf_compare_output]

'''
'''
vcf.infos => ordered dict of {'TAG' : vcf.parser.Info}
vcf will automatically handle the type=Interger, etc. field for us.
Note that ALT might be multiple alleles
'''
import pandas as pd
import operator
from functools import partial

'''
#TODO: Replace NaN with None ala:
df.where(pd.notnull(df), None)
df.fillna("NULL")
'''

def match_records(left, right, tag, threshold):
    l, r = flatten_vcf(left), flatten_vcf(right)
    assert tag in l and tag in r, "Tag not found in flattened record {0}".format(str(flatten_vcf))
    if not threshold:
        return l[tag] == r[tag]
    else:
        return abs(l[tag] - r[tag]) < threshold

def vcall(records):
    return [rec for rec in records if rec.REF != rec.INFO['CB']]

def diff(left, right, tag, threshold):
    '''
    :param list left: vcf Records
    :param list right: vcf Records
    :param str tag: i.e. CB, ALT
    :param int thershold: the amount of difference necessary
    :return list of tuples from left, right where the values differed
    '''
    #return [(l, r) for l, r in zip(left, right) if not match_records(l, r, tag, threshold)]
    return [l for l, r in zip(left, right) if not match_records(l, r, tag, threshold)]

def validate_vcf(filename, tags):
    #TODO: Unfortunately have to validate types here.
    # Also Check that tags are in vcf (using flattenvcf)
    pass

def compare_value(comparator, val, rec, tag):
    flat_vcf = flatten_vcf(rec)
    assert tag in flat_vcf, "Tag not found in flattened record {0}".format(str(flatten_vcf))
    return comparator(flat_vcf[tag], val)

def not_members(vcf_list, tag, collection):
    return [rec for rec in vcf_list if flatten_vcf(rec)[tag] not in collection]

ambiguous = partial(not_members, tag='CB', collection=['A', 'C', 'G', 'T'])
exists = partial(not_members, collection= [None, [None], '-'])


def filter(vcf_list, tag, opfunc, value):
    '''
    i.e.  CBD > 12, etc.
    :param list vcf_list: list of vcf.model._Record objects as returned from vcf.Reader
    :param str operator: a function contained in the `operator` python stdlib module
    :param str tag: a valid field in the vcf file (i.e. ALT, CBD)
    :param object value: str, int/float, or list. the value to run operator against
    :return a subset of a vcf record where the operator evaluates to true
    '''
    #assert ( op in ['exists', 'ambiguous']  != bool(value)), "filter should not be called with operator 'exisits' OR a value."
    #ALTs are stored as lists by vcf
    if tag == 'ALT' and type(value) is str:
        value = [base.strip()  for base in value.split(',')]
    if type(opfunc) is str:
        opfunc = getattr(operator, opfunc)
    check_value = partial(compare_value, opfunc, value)
    return [rec for rec in vcf_list if check_value(rec, tag)]


#Can use the dtype of a column
''' file.samples()[0]  -> returns dengue.bam etc.'''

def flatten_list(A):
  return A[0] if type(A) == list else A

def flatten_vcf(record):
      _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT']
      fields = [getattr(record, _id) for _id in _ids]
     # converters = [str, str, str,  int, flatten_list]
     # fields = (convert(getattr(record, _id)) for convert, _id in zip(converters, _ids))
      d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
      d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
      return d

def vcf_file_to_df(vcf_records):
    return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)

def stat(records):
   return vcf_file_to_df(records).describe()

def statdiff(recs_a, recs_b):
   df1, df2 = vcf_file_to_df(recs_a), vcf_file_to_df(recs_b)
   return df1.describe() - df2.describe()

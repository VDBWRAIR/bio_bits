'''
Usage:
    vcf_compare <vcfiles> [--out-directory=<outdir>] 
                          [--threshold=THRESH]
                          [--fields=FIELDS]...

Options:
    -o DIR --output=DIR  The output directory [default: vcf_compare_output]
    -t THRESHOLD --threshold=THRESHOLD  Numerical threshold to determine difference
    --


vcf.infos => ordered dict of {'TAG' : vcf.parser.Info}
vcf will automatically handle the type=Interger, etc. field for us. 
Note that ALT might be multiple alleles (sadface)      


'''
import vcf
import pandas as 
import itertools

'''
Flatten multiple fastq files into one iterator.
'''
def multi_fastq_iterator(fastq_handles):
    
    files = map(lambda a: SeqIO.parse(open(a), 'fastq'), fastq_handles)
    #? return itertools.chain(files)
    return record for _fastq in files for record in _fastq
    

''' file.samples()[0]  -> returns dengue.bam etc.''' 
def SmartVCF(vcf.parser.Reader):
    def __init__(self, *args, **kwargs):
        return super(vcf.parser.Reader, self).__init__(*args, **kwargs)

def flatten_list(A): 
  return A[0] if type(A) == list else A

def flatten_vcf(record):
      _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT']
      fields = getattr(record, _id) for _id in _ids
     # converters = [str, str, str,  int, flatten_list]
     # fields = (convert(getattr(record, _id)) for convert, _id in zip(converters, _ids))
      d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
      d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
      return d

def vcf_file_to_df(filename):
   vcf_records = vcf.Reader(open(filename))
   return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)

def get_statistics_diff(file_a, file_b):
   df1, df2 = vcf_file_to_df(file_a, file_b)
   return df1.describe() - df2.describe()
'''
i.e.
sum(1 for base in sub_df.ALT if base)
sum(1 for base1, base2 in zip(sub_df.CB, norms.CB) if base1 != base2) 
sum(1 for base1, base2 in zip(sub_df.HPOLY, norms.HPOLY) if base1 != base2)
'''
def get_num_different(df1, df2, column):
    pass


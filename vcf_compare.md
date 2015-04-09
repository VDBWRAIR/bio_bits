



Script will expect The following arguments:
1. two .vcf files
2.  those fields to check for differences on (i.e. Called base, etc.). Optional, check all fields by default
3. A numerical (percent-based?) thershold for deciding if two attributes are different. Optional, 0.01 by default
   i.e. CBD at position 0 is:
      file1: 12
      file2: 14
if threshold is 2 or less, than the record at position 0 will be included under the `CBD` heading in the `diffs.txt` file.



I propose the ouoptut will be in the form of the following files:
`sumarry.txt`

Statistical summary of the differences between the two files, like average and total call depth difference, etc.

`diffs.txt`
A field by field listing of how many different records appeared for each field, followed by the differing vcf records (as seen in the original files) listed with their position.
i.e., file1 has CBD=99 at position 0, while file2 has CBD=150 at position 0, similar at position 9
```
There were 2 records with differing called base depth with threshold 10:
0> VCF record from file 1
0< VCF record from file 2
9< VCF record from file 1
9> VCF record from file 2
```
`headers.txt`
The ouptut of `diffs.txt` without the actual records listed. i.e.:
```
There were 2 records with differing called base depth with threshold 10:
There were 4 records with differing called bases.
```
`diffbases.diff`
The output of the existing `vcf_diff.py` script (showing where Called bases were different from references) piped to `diff`.  Or, these records listed side-by-side, as part of `diffs.txt`.

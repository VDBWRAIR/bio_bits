import pandas as pd
import numpy as np
from StringIO import StringIO


# Simple Nucleotide Substitution matrix
# Copy/pasted table from 
# http://www.jalview.org/help/html/calculations/scorematrices.html#simplenucleotide
sns = StringIO('''
,A,C,G,I,N,R,T,U,X,Y
A,10,-8,-8,1,1,1,-8,-8,1,-8
C,-8,10,-8,1,1,-8,-8,-8,1,1
G,-8,-8,10,1,1,1,-8,-8,1,-8
I,1,1,1,10,1,0,1,1,0,0
N,1,1,1,1,10,1,1,1,1,1
R,1,-8,1,0,1,10,-8,-8,0,-8
T,-8,-8,-8,1,1,-8,10,10,1,1
U,-8,-8,-8,1,1,-8,10,10,1,1
X,1,1,1,0,1,0,1,1,10,0
Y,-8,1,-8,0,1,-8,1,1,0,10
''')
# Converted to pandas dataframe indexed in both directions
sns_matrix = pd.read_csv(sns, index_col=0)

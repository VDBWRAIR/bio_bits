'''
Usage: make_pca.py <fasta> [--map <mapfile>] [--outdir <DIR>] [--coord <coordfile>]

Options:
    --outdir=<DIR>,-o=<DIR>   Directory to put html file in. [Default: emperor]
    --map=<mapfile>,-m=<mapfile>   TSV file which maps FASTA IDs to metadata. If not supplied one is generated using the FASTA IDs only.
    --coord=<coordfile>,-c=<coordfile>   Coordinate file including distance matrix, defined by Qiime pipeline
'''
#options in docopt are special, and need = if using them
''' also create static png files.
allow "color-by" parameters. '''
from skbio import Alignment
from skbio.stats.ordination import PCoA
from itertools import groupby
from docopt import docopt
from schema import Schema, Use, Optional
import os
import sh
import re
import emperor #!not used, but emperor must be installed!

run_emperor = sh.Command('make_emperor.py')

def make_coordinates(fasta_filename):
    alignment = Alignment.read(fasta_filename)
    distance_matrix = alignment.distances()
    pcoa = PCoA(distance_matrix)
    scores = pcoa.scores()
    return scores

def write_coordiates(fasta_filename):
    outname = '%s.coord' % fasta_filename
    assert not os.path.exists(outname), "Coordinate file %s exists! Please remove or run again with --coord parameter." % outname
    make_coordinates(fasta_filename).write(outname)
    return outname

def make_emperor(fasta_fn, outdir, mapfile, coordfile):
    # do try/except for sh call to make_emperor
    mapfile = mapfile or make_simple_mapping(fasta_fn)
    coordinate_file = coordfile or write_coordiates(fasta_fn)
    return run_emperor(i=coordinate_file, m=mapfile, o=outdir)

def make_simple_mapping(fasta_fn):
    ids = map(lambda x: x[1:], filter(lambda x: x.startswith('>'), open(fasta_fn)))
    header = '#SampleID\n'
    mapfile_fn = '%s.map' % fasta_fn
    assert not os.path.exists(mapfile_fn), "Mapping file %s exists! Please remove, or run again with --map parameter." % mapfile_fn
    with open(mapfile_fn, 'w') as mapfile:
        mapfile.write(header)
        mapfile.writelines(ids)
    return mapfile_fn

#NOTE: Currently unused
'''
def make_undescore_metadata_mapping(fasta):
    reg = re.compile(r'^[^_]+_([^_]+)_')
    with open('%s.map' % fasta, 'w') as mapfile:
        ids = map(X[1:], filter(X[0] == '>', open(fasta)))
        #groups = groupby(ids, lambda x: x[:x.find('_')])
        groups = groupby(ids, lambda x: reg.match(x).groups()[0])
        header = '#Group\tSampleID\n'
        mapfile.write(header)
        for k, group in groups:
            mapfile.writelines(map(('%s\t' % k + '{0}').format, group))
    return mapfile.name
'''
def main():
    scheme = Schema(
        { '<fasta>' : os.path.isfile,
         Optional('--map') : Use(lambda x: x or None),
         Optional('--coord') : Use(lambda x: x or None),
         '--outdir' : str
         })
    raw_args = docopt(__doc__, version='Version 1.0')
    args = scheme.validate(raw_args)
    make_emperor( args['<fasta>'], args['--outdir'], args['--map'], args['--coord'])

if __name__ == '__main__': main()

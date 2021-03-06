'''
Usage:
    plot_muts.py --query <query> --refs <refs> [--out <outfile>] [--html] [--cluster]

Options:
    --cluster            Evaluate by comparing to the first two sequences in <refs> rather than by time
    --html               Produce html output in addition to a static image
    --refs,-r=<refs>     Fasta file, sequence with earliest year is base reference
    --query,-q=<query>   Query sequences
    --out,-o=<outfile>   Figure saved here

Help:
    All sequences must be the same length.
'''
from __future__ import print_function
import numpy as np
from functools import partial
import operator
import os, sys, re
from Bio import SeqIO
import matplotlib.pyplot as plt
import docopt, schema
from operator import itemgetter as get
import csv
from dateutil import parser
import datetime
from time import mktime
from funcy import compose
from funcy.py2 import map, zip
try:
    #below import is necessary for some reason
    from scipy.stats import poisson
    import scipy
    DISTRIBUTION = scipy.stats.poisson
except ImportError:
    DISTRIBUTION = None
years = r'190\\d|19[1-9]\\d|200\\d|201[0-5]' #1900-2015
year_regex = re.compile(years)
hamming = compose(sum, partial(map, operator.ne))
timestamp = lambda x: mktime(x.timetuple())
legend = {"queries": 'r', "references": 'b', "interval": 'g'}
#def pdist(s1, s2):
#    assert len(s1) == len(s2), "All sequences must be the same length! %s %s" % (s1, s2)
#    return hamming(s1, s2)/float(len(s1))
class InvalidFastaIdentifier(Exception): pass
def extract_date(fasta_id):
    ''' fasta id '''
    _e = InvalidFastaIdentifier("Could retrieve date from {0}".format(fasta_id))
    if '____' not in fasta_id:
        raise _e
    s = fasta_id.split('____')[-1]
    try:
        dt = parser.parse(s.replace('_','/'))
        return dt
    except Exception as e:
        print("Error parsing {0}".format(s))
        raise _e

def get_seqs_and_dates(fn):
    fasta = SeqIO.parse(fn, format="fasta")
    info = [(str(seq.seq), seq.id, seq.description) for seq in fasta]
    seqs, ids, descriptions = zip(*info)
    dates = map(extract_date, ids)
    return seqs, dates, ids

def process(refs_fn, query_fn, save_path=None, html=True):
    ref_seqs, ref_dates, ref_names = zip(*sorted(zip(*get_seqs_and_dates(refs_fn)), key=get(1)))
    #assert len(ref_seqs) > 1, "Need more than 1 reference sequence"
    ref_seqs = map(str.upper, ref_seqs)
    super_ref_seq, super_ref_date, super_ref_name = ref_seqs[0], ref_dates[0], ref_names[0]
    get_mutations = partial(hamming, super_ref_seq)
    def get_relative_info(seqs, dates, names):
        muts = map(get_mutations, seqs)
        dists = [(yr - super_ref_date).days for yr in dates]
        return muts, dists, names
    ref_muts, ref_dists, ref_names =  get_relative_info(ref_seqs, ref_dates, ref_names)
    query_muts, query_dists, query_names = get_relative_info(*get_seqs_and_dates(query_fn))
    do_plot(ref_dists, ref_muts, ref_names, query_dists, query_muts, query_names, save_path, html)
    #map(compose(print, '{0}\t{1}'.format ), ref_dists, ref_muts)

def do_plot(x1, y1, ref_names, x2, y2, query_names, save_path=None, html=True, \
            title='Mutations over time (days)', x_axis_label='time since base reference', y_axis_label='p-distance'):
    '''
    :param iterable x1: reference dates distances
    :param iterable y1: reference p-distances
    :param iterable x2: query dates diferences
    :param iterable y2: query p-distances
    :param str save_path: path to save image or None to open GTK if installed
    '''
    assert len(x1) > 0, "No reference dates to use"
    assert len(y2) > 0, "No reference p-distances to use"
    assert len(x2) > 0, "No query dates to use"
    assert len(y2) > 0, "No query p-distances to use"
    assert len(ref_names) == len(x1) and len(query_names) == len(x2)
    fig = plt.figure()
    ax = plt.subplot(111)
#    from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
#    years = YearLocator()   # every year
#    months = MonthLocator()  # every month
#    yearsFmt = DateFormatter('%Y')
#    ax.xaxis.set_major_locator(years)
#    ax.xaxis.set_major_formatter(yearsFmt)
#    ax.xaxis.set_minor_locator(months)
    max_x = max(max(x1), max(x2))
    #legend_info = [mpatches.Patch(label=n, color=c) for n, c in legend.items()]
    """ http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot"""

    ref_info = zip(ref_names, x1, y1)
    query_info = zip(query_names, x2, y2)
    all_info = sorted(ref_info + query_info, key=lambda x: x[2], reverse=True)

    if save_path:
        fh = open(save_path+'.csv', 'wt')
    else:
        fh = sys.stdout
    fh.write('name,dates,p-dist\n')
    outcsv = csv.writer(fh)
    map(outcsv.writerow, all_info)
    if html:
        assert sys.version[:3] != '2.6', "Requires python 2.7 or higher to run bokeh."
        import bokeh.models as bkm
        import bokeh.plotting as bkp
        bokeh_tools = [bkm.WheelZoomTool(), bkm.PanTool(), bkm.BoxZoomTool(),
             bkm.PreviewSaveTool(), bkm.ResetTool(), bkm.BoxSelectTool(),
             bkm.ResizeTool()]
        ref_names = map('R: {0}'.format, ref_names)
        query_names = map('Q: {0}'.format, query_names)
        hover = bkm.HoverTool(tooltips=[("id", "@ids"),]) #   ("(days,muts)", "($x, $y)"),
        source1 = bkm.ColumnDataSource(data=dict(x=x1, y=y1,  ids=ref_names))
        source2 = bkm.ColumnDataSource(data=dict(x2=x2, y2=y2, ids=query_names))
        p = bkp.figure(plot_width=400, plot_height=400, tools=[hover]+bokeh_tools, title=title)
        p.circle('x', 'y', source=source1, line_color='gray', legend='reference')
        p.square('x2', 'y2', source=source2, fill_color='red', legend='query')
        if save_path:
            bkp.output_file(save_path + '.html')
            bkp.save(p)
        else: bkp.show(p)
    else:
        plot_muts(ax, x1, y1, plotkwargs=dict(label='references (blue)', color=legend['references'], marker='s'), polyfit=True, max_x=max_x, dist=None)
        plot_muts(ax, x2, y2, plotkwargs=dict(label='queries (red)', color=legend['queries']), dist=None)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(handles=legend_info, loc='center left', bbox_to_anchor=(1, 0.5))
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), framealpha=0)
        ax.legend(framealpha=0)
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()

def plot_muts(ax, x, y, dist=DISTRIBUTION, polyfit=False, max_x=None, plotkwargs=dict(marker='o')):
    '''
    Plot x and y

    if norm distribution, probably have to scale (via passing loc= and scale=)
    problem was didn't account for +b
    '''
    if max_x and isinstance(max_x, datetime.datetime):
        max_x = timestamp(max_x)
    retval = ax.scatter(x, y, **plotkwargs)#color=color, label=label, marker=marker)
    if polyfit:
        ''' this forces a polyfit with y-intercept at zero, necessary because
        we necessarily start with 0 mutations from the query at date 0.'''
        x = np.array(x)[:,np.newaxis]
        m, _, _, _ = np.linalg.lstsq(x, y)
        x, y = np.linspace(0,max_x,100), m*np.linspace(0,max_x,100)
        #ax.plot(x, y, color='y', label='Best Fit', linewidth=2)
        ax.plot(x, y, color='y', linewidth=2)
    if dist:
        """see  http://stackoverflow.com/a/14814711/3757222"""
        R = dist.interval(0.95, y)
        interval_left, interval_right = R
        interval_color = legend['interval']
        ax.plot(x, interval_left, color=interval_color)
        ax.plot(x, interval_right,color=interval_color)
    return retval

#def test_more():
#    refs = range(25), range(25)
#    queries = [1, 5, 20, 10], [2, 20, 40, 10]
#    do_plot(refs[0], refs[1], queries[0], queries[1], None)
#
#def test_plot():
#    ''' can verify this works by using scipy.stats.norm.interval instead'''
#    default_x = range(25)
#    default_y = range(0, 50, 2)
#    plot_muts(default_x, default_y, 'r', True, scipy.stats.poisson, max_x=max(default_x))
#    plt.show()

def extract_info(fasta):
    objs = SeqIO.parse(fasta, format="fasta")
    info = [(str(seq.seq), seq.id) for seq in objs]
    seqs, ids = zip(*info)
    return seqs, ids

def get_clusters(refs, queries):
    all_ref_seqs, all_ref_ids = extract_info(refs)
    ref1, ref2 = all_ref_seqs[:2]
    dists1, dists2 = partial(hamming, ref1), partial(hamming, ref2)
    ref_seqs, ref_ids = all_ref_seqs[2:], all_ref_ids[2:]
    # There may only be 2 references to compare
    # so we check here if that is the case and set distances to 0 if so
    if ref_seqs:
        ref_dists1, ref_dists2 = map(dists1, ref_seqs), map(dists2, ref_seqs)
    else:
        ref_dists1 = ref_dists2 = [0,0]
        ref_ids = ['','']

    query_seqs, query_ids = extract_info(queries)
    query_dists1, query_dists2 = map(dists1, query_seqs), map(dists2, query_seqs)

    return ref_dists1, ref_dists2, ref_ids, query_dists1, query_dists2, query_ids, all_ref_ids[0], all_ref_ids[1]

def process_cluster(refs, queries, save_path=None, html=False):
    ref_dists1, ref_dists2, ref_ids, query_dists1, query_dists2, query_ids, ref1_id, ref2_id = get_clusters(refs, queries)
    do_plot(ref_dists1, ref_dists2, ref_ids, query_dists1, query_dists2, query_ids, \
            save_path, html, title='test', x_axis_label=ref1_id, y_axis_label=ref2_id)

def main():
    #if sys.argv[1] == 'test': test_more()
    scheme = schema.Schema(
        { '--query' : os.path.exists,
          '--refs' : os.path.exists,
         schema.Optional('--out') : lambda x: True,
         schema.Optional('--html') : lambda x: True,
         schema.Optional('--cluster') : lambda x: True
        # schema.Or(lambda x: x is None,  #check file can be created
        #                                      lambda x: os.access(os.path.dirname(x), os.W_OK))
         })
    args = docopt.docopt(__doc__, version='Version 1.0')
    scheme.validate(args)
    queries, refs, out, do_cluster = args['--query'], args['--refs'], args['--out'], args['--cluster']
    if do_cluster:
        process_cluster(refs, queries, out, args['--html'])
    else:
        process(refs, queries, out, args['--html'])

if __name__ == '__main__': main()

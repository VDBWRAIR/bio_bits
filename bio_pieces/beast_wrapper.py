import sys
import re

import sh

import beast_est_time

def get_xmlpath_from_argv(argv):
    r'''
    Get the beast input xml file from argv list

    >>> argv = ['beast', 'foo.xml', '-beagle_GPU']
    >>> get_xmlpath_from_argv(argv)
    'foo.xml'
    >>> argv[1] = '/path/to/foo.xml'
    >>> get_xmlpath_from_argv(argv)
    '/path/to/foo.xml'
    >>> get_xmlpath_from_argv(argv[:1]+argv[-1:])
    ValueError: Missing input xml file path
    '''
    try:
        xmlpath = filter(lambda x: x.endswith('.xml'), argv)[0]
        return xmlpath
    except IndexError:
        raise ValueError('Missing input xml file path')

def get_chainlength_from_xml(xml_fh):
    r'''
    Get chainLength from opened xml file handle

    >>> from StringIO import StringIO
    >>> xml = StringIO('foo bar baz\n"<\n chainLength="10000" bar\n baz')
    >>> get_chainlength_from_xml(xml)
    '10000'
    >>> get_chainlength_from_xml(StringIO('foo\nbar baz'))
    ValueError
    '''
    chainlength = None
    p = re.compile('chainLength="(\d+)"')
    for line in xml_fh:
        m = p.search(line)
        if m:
            chainlength = m.group(1)
    if chainlength is None:
        raise ValueError('Could not get chainLength from {0}'.format(xml_fh))
    return chainlength

def run_beast(*args, **kwargs):
    '''
    Simply run beast with same args supplied to beast_wrapper but add sec_to_time
    column for remaining time left
    '''
    xmlpath = get_xmlpath_from_argv(args)
    chainlength = get_chainlength_from_xml(open(xmlpath))
    for line in sh.beast(*args, _iter=True):
        if 'hours/' in line:
            sec = beast_est_time.hours_states_to_sec(line, chainlength)
            pretty_time = beast_est_time.sec_to_time(sec)
            sys.stdout.write(line.rstrip() + '\t' +  pretty_time + '\n')
        else:
            sys.stdout.write(line)

def main():
    run_beast(*sys.argv[1:])

if __name__ == '__main__':
    main()

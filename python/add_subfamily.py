#! /usr/bin/env python

import sys
import re
import utils
from collections import defaultdict

def main(args):
    ltrtypes = defaultdict(list)
    gtf = utils.tab_line_gen(args.gtf)
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        if attrd['repType']=='ltr':
            ltrtypes[attrd['locus']].append(attrd['repName'])
    
    titer = utils.tab_line_gen(args.infile)
    header = titer.next()
    print >>args.outfile, '\t'.join(header + ['subfamily'])

    for row in titer:
        if row[0] in ltrtypes:
            subfam  = ','.join(sorted(set(ltrtypes[row[0]])))
        else:
            subfam = ''
        
        print >>sys.stdout, '\t'.join(row + [subfam])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Add locus field to each annotation.')
    parser.add_argument('--gtf', type=argparse.FileType('rU'), required=True)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())

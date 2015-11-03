#! /usr/bin/env python

import utils
from collections import defaultdict
import re

def main(parser):
    args = parser.parse_args()
    
    rows = []
    base_columns = ['chrom','source','feature','start','end','score','strand','frame']
    attr_columns = set()
    
    # Parse rows
    gtflines = utils.tab_line_gen(args.infile)
    for l in utils.sort_gtf(gtflines):
        rowd = dict(zip(base_columns, l[:8]))
        for k,v in re.findall('(\S+)\s+"([\s\S]+?)";',l[8]):
            attr_columns.add(k)
            rowd[k] = v
        rows.append(rowd)
    
    # Set column headers
    columns = base_columns + list(attr_columns)
    if args.keycol in columns:
        columns.remove(args.keycol)
        columns = [args.keycol] + columns
    
    # Print the table
    print >>args.outfile, '\t'.join(columns)
    for rowd in rows:
        print >>args.outfile, '\t'.join([rowd[c] if c in rowd else '' for c in columns])

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Output table from GTF')
    parser.add_argument('--keycol', default='locus')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser)

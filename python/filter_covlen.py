#! /usr/bin/env python

import sys
import os
from collections import defaultdict, Counter
import re
import utils

def main(args):
    discard = set()
    discard_fh = open(args.discard, 'w')

    lines = utils.tab_line_gen(args.infile)
    for l in lines:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        if l[1] == 'merged':
            cov = int(attrd['cov'])
            if cov < args.threshold:
                discard.add(attrd['name'])
                print >>discard_fh, '\t'.join(l)
            else:
                print >>args.outfile, '\t'.join(l)            
        else:
            if attrd['gene_id'] in discard:
                print >>discard_fh, '\t'.join(l)
            else:
                print >>args.outfile, '\t'.join(l)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--threshold', type=int, default=500)
    parser.add_argument('--discard', default=os.devnull)    
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())
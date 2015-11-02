#! /usr/bin/env python

import sys
import re
import utils
from collections import defaultdict

def main(args):
    if args.mapping:
        tdict = {}
        table = utils.tab_line_gen(args.mapping)
        header = table.next()
        for t in table:
            d = dict(zip(header,t))
            tdict[d[args.key]] = d
    else:
        # Lookup of any key in tdict returns a dictionary where the value of 'name' is ''
        tdict = defaultdict(lambda:{'name':''})

    gtf = utils.tab_line_gen(args.infile)
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        assert 'locus' not in attrd    
        if g[1]=='merged':
            newname = tdict[attrd['name']][args.value]
            if newname == '':
                n1,n2 = attrd['name'].split('_')
                newname = '%s_%s_%s' % (n1, g[0].strip('chr'), n2)
            g[8] += ' locus "%s";' % newname
            print >>args.outfile, '\t'.join(g)
        else:
            newname = tdict[attrd['gene_id']][args.value]
            if newname == '':
                n1,n2 = attrd['gene_id'].split('_')
                newname = '%s_%s_%s' % (n1, g[0].strip('chr'), n2)
            g[8] += ' locus "%s";' % newname
            print >>args.outfile, '\t'.join(g)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Add locus field to each annotation.')
    parser.add_argument('--mapping', type=argparse.FileType('rU'), help='File mapping annotation IDs to names')
    parser.add_argument('--key', default='id', help='Column in mapping file with annotation ID')
    parser.add_argument('--value', default='name', help='Column in mapping file with replacement value')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())

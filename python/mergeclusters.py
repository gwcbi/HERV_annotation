#! /usr/bin/env python
import sys
from collections import defaultdict, Counter
import re
import utils

def main(parser):
    args = parser.parse_args()
    lines = utils.tab_line_gen(args.infile)
    clustered = utils.cluster_gtf(lines)
    reptypes = utils.by_attribute(clustered, 'repType')

    catcount = Counter()
    newlines = []
    for cnum,c in clustered.iteritems():
        c.sort(key=lambda x:int(x[3]))
        #
        if reptypes[cnum] == ['ltr','internal','ltr']:
            category = 'canonical'
        elif reptypes[cnum] == ['ltr','internal'] or reptypes[cnum] == ['internal','ltr']:
            category = 'oneside'
        elif reptypes[cnum] == ['internal']:
            category = 'soloint'
        elif reptypes[cnum] == ['ltr']:
            category = 'sololtr'
        else:
            category = 'unusual'
        catcount[category] += 1
        #
        cstart = min(int(l[3]) for l in c)
        cend   = max(int(l[4]) for l in c)
        strands = set(l[6] for l in c)
        if len(strands)==1: cstrand = strands.pop()
        else: cstrand = '.'
        cattr = 'name "%s_%04d"; category "%s";' % (args.prefix, int(cnum), category)
        cann = [c[0][0], 'merged', 'gene', str(cstart), str(cend), '.', cstrand, '.', cattr]
        newlines.append(cann)

        for l in c:
            l[1] = category        
            attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
            if 'gene_id' in attr: del attr['gene_id']
            if 'transcript_id' in attr: del attr['transcript_id']
            l[8] = 'gene_id "%s_%04d"; transcript_id "cluster%04d"; ' % (args.prefix, int(cnum), int(cnum))
            l[8] = l[8] + ' '.join('%s "%s";' % (k,v) for k,v in attr.iteritems())
            newlines.append(l[:-1])

    for l in utils.sort_gtf(newlines):
        print >>args.outfile, '\t'.join(l)       

    for cat in ['canonical','oneside','soloint','sololtr','unusual']:
        print >>sys.stderr, '%s:     %d' % (cat, catcount[cat])

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='process clustered annotations')
    parser.add_argument('--prefix', default='cluster')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser)

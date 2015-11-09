#! /usr/bin/env python

import sys
import re
import utils
from collections import defaultdict

def main(args):
    names = args.names.split(',')
    category = args.category
    gtf = utils.tab_line_gen(args.infile)

    # Skip lines not being merged
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        if 'name' in attrd and attrd['name'] in names:
            break
        else:
            print >>args.outfile, '\t'.join(g)

    mergelines = [g]
    sublines   = []
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        if g[1]=='merged' and attrd['name'] in names:
            mergelines.append(g)
        elif 'gene_id' in attrd and attrd['gene_id'] in names:
            sublines.append(g)
        else:
            break

    # Merge the merge lines    
    new_merge = mergelines[0][:]
    nm_spos = min(int(m[3]) for m in mergelines)
    nm_epos = max(int(m[4]) for m in mergelines)
    nm_strand = set(m[6] for m in mergelines)
    
    all_attrd = [dict(re.findall('(\S+)\s+"([\s\S]+?)";', m[8])) for m in mergelines]
    nm_name = all_attrd[0]['name']
    nm_category = category if category is not None else all_attrd[0]['category']
    nm_nfeats = sum(int(d['nfeats']) for d in all_attrd)
    nm_length = (nm_epos - nm_spos)
    nm_cov = sum(int(d['cov']) for d in all_attrd)
    
    new_merge[3] = str(nm_spos)
    new_merge[4] = str(nm_epos)
    new_merge[6] = nm_strand.pop() if len(nm_strand) == 1 else '.'
    new_merge[8] = 'name "%s"; category "%s"; nfeats "%d"; length "%d"; cov "%d";' % (nm_name, nm_category, nm_nfeats, nm_length, nm_cov)
    print  >>args.outfile, '\t'.join(new_merge)

    for l in sublines:
        l[1] = category if category is not None else all_attrd[0]['category']
        l[8] = re.sub('gene_id "\S+";', 'gene_id "%s";' % nm_name , l[8])
        l[8] = re.sub('transcript_id "\S+";', 'transcript_id "%s";' % nm_name , l[8])
        print >>args.outfile, '\t'.join(l)

    # Resume printing the file
    print >>args.outfile, '\t'.join(g) 
    for g in gtf:
        print >>args.outfile, '\t'.join(g)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Add locus field to each annotation.')
    parser.add_argument('--names', help='Comma seperated list of names to merge')
    parser.add_argument('--category', help='Category for merged annotation')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())

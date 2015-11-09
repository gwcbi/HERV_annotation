#! /usr/bin/env python

import sys
import re
import utils
from collections import defaultdict

# Make new merge line
def make_mergeline(clust, name, category):
    nm_spos     = min(int(l[3]) for l in clust)
    nm_epos = max(int(l[4]) for l in clust)
    nm_strand = set(l[6] for l in clust)
    nm_strand = nm_strand.pop() if len(nm_strand) == 1 else '.'
    nm_cov = utils.covered_len(clust)
    nm_attr = 'name "%s"; category "%s"; nfeats "%d"; length "%d"; cov "%d";' % (name, category, len(clust), (nm_epos-nm_spos), nm_cov)
    nm = [clust[0][0], 'merged', 'gene', str(nm_spos), str(nm_epos), '.', nm_strand, '.', nm_attr]
    return nm


def main(args):
    # Load GTF file
    gtf = utils.tab_line_gen(args.infile)

    # Skip lines not being merged
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        if 'name' in attrd and attrd['name'] == args.name:
            break
        else:
            print >>args.outfile, '\t'.join(g)

    # Get all lines for locus to be split
    mergelines = [g]
    sublines   = []
    for g in gtf:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        if g[1]=='merged' and attrd['name'] == args.name:
            mergelines.append(g)
        elif 'gene_id' in attrd and attrd['gene_id'] == args.name:
            sublines.append(g)
        else:
            break

    # Set the new category and names
    if args.category:
        category1, category2 = args.category.split(',')
    else:
        attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',g[8]))
        category1 = attrd['category']
        category2 = attrd['category']

    if args.newname:
        newname1, newname2 = args.newname.split(',')
    else:
        newname1 = '%s.1' % args.name
        newname2 = '%s.2' % args.name
    
    # Print out the edited lines
    sub1 = sublines[:args.split]
    print >>args.outfile, '\t'.join(make_mergeline(sub1, newname1, category1))
    for l in sub1:
        l[1] = category1
        l[8] = re.sub('gene_id "\S+";', 'gene_id "%s";' % newname1, l[8])
        l[8] = re.sub('transcript_id "\S+";', 'transcript_id "%s";' % newname1, l[8])
        print >>args.outfile, '\t'.join(l)

    sub2 = sublines[args.split:]
    print >>args.outfile, '\t'.join(make_mergeline(sub2, newname2, category2))
    for l in sub2:
        l[1] = category2
        l[8] = re.sub('gene_id "\S+";', 'gene_id "%s";' % newname2, l[8])
        l[8] = re.sub('transcript_id "\S+";', 'transcript_id "%s";' % newname2, l[8])
        print >>args.outfile, '\t'.join(l)

    # Resume printing the file
    print >>args.outfile, '\t'.join(g) 
    for g in gtf:
        print >>args.outfile, '\t'.join(g)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Add locus field to each annotation.')
    parser.add_argument('--name', help='Name of locus to split')
    parser.add_argument('--split', type=int, help='Index for split')
    parser.add_argument('--newname', help='New names for loci')
    parser.add_argument('--category', help='Categories for new loci')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())

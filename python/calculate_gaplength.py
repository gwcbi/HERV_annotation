#! /usr/bin/env python
import re
import utils
import sys

def main(parser):
    args = parser.parse_args()
    lines = utils.tab_line_gen(args.infile)
    bystrand = {'+':[], '-':[]}
    for l in lines:
        bystrand[l[6]].append(l)

    bystrand['+'] = list(utils.sort_gtf(bystrand['+']))
    bystrand['-'] = list(utils.sort_gtf(bystrand['-']))
    
    grouped = {'+':[], '-':[]}
    for strand in ['+','-']:
        score = None
        chrom = None
        tmp = []    
        for l in bystrand[strand]:
            if score is not None:
                if l[5] != score or l[0] != chrom:
                    grouped[strand].append(tmp)
                    tmp = []
            tmp.append(l)
            score = l[5]
            chrom = l[0]

    gaplens = []
    merged = []
    for g in grouped['+'] + grouped['-']:
        if len(g) == 1:
            merged.append(g[0])
        else:
            s = ''
            for i in range(len(g)-1):
                gaplen = int(g[i+1][3]) - int(g[i][4])
                s += '%s:%s-%s(%s)' % (g[i][0],g[i][3],g[i][4],g[i][6])
                s += ' --- %d --- ' % gaplen
                gaplens.append(gaplen)
                
            s += '%s:%s-%s(%s)' % (g[-1][0],g[-1][3],g[-1][4],g[-1][6])
            print >>sys.stderr, s
            # spos = min(int(l[3]) for l in g)
            # epos = max(int(l[4]) for l in g)
            # attrs = [dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8])) for l in g] 
            # newline = [g[0][0], 'joined', 'exon', str(spos), str(epos), g[0][5], g[0][6], '.']
            # newattr = {'joined': ','.join(a['id'] for a in attrs),
            #            'repType': attrs[0]['repType'],
            #            }
            # newline.append(' '.join('%s "%s";' % (k,v) for k,v in newattr.iteritems()))
            # merged.append(newline)

    if gaplens:
        print >>sys.stderr, 'min gap length:  %d' % min(gaplens)
        print >>sys.stderr, 'mean gap length: %d' % (float(sum(gaplens)) / len(gaplens))
        print >>sys.stderr, 'max gap length:  %d' % max(gaplens)
    else:
        print >>sys.stderr, 'No gaps found'

    print >>args.outfile, '%d' % max(gaplens)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Calculates length of gaps between annotations that are part of the same hit')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser)

#! /usr/bin/env python

import sys
import re
import utils
from collections import defaultdict

cyto = dict([l.strip('\n').split('\t') for l in open('tmp/cyto_name_map.txt','rU')])
lines = list(utils.tab_line_gen(open('filtered.hg19.gtf','rU')))

### Seperate merged and unmerged lines 
### Group unmerged lines according to gene
merged = []
unmerged = defaultdict(list)
for l in lines:
    attrs = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
    if l[1]=='merged':
        merged.append((attrs['name'],l))
    else:
        attrs = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        unmerged[attrs['gene_id']].append(l)

sdata = defaultdict(dict)
for mname,ml in merged:
    attrs = dict(re.findall('(\S+)\s+"([\s\S]+?)";',ml[8]))
    sdata[mname]['category'] = attrs['category']
    sdata[mname]['band'] = cyto[mname]
    sdata[mname]['id'] = mname
    sdata[mname]['start'] = int(ml[3])
    sdata[mname]['chrom'] = ml[0]
    sdata[mname]['alias'] = ''

bandorder = []
byband = defaultdict(list)
for mname in sorted(sdata.keys()):
    d = sdata[mname]
    if d['category'] != 'sololtr':
        byband[d['band']].append(d)
        if d['band'] not in bandorder: bandorder.append(d['band'])
    elif d['alias']:
        byband[d['band']].append(d)
        if d['band'] not in bandorder: bandorder.append(d['band'])        

someletters = 'abcdefghijklmnopqrstuvwxyz'
for band in bandorder:
    dlist = byband[band]
    dlist.sort(key=lambda x:x['start'])
    if len(dlist) == 1:
        d = dlist[0]
        sdata[d['id']]['name'] = 'HERV9_%s' % band
        if d['alias'] != '' and d['alias'] != band:
            print >>sys.stderr, 'Mismatch: %s %s' % (d['alias'],band)
            d['alias'] += '*'
    else:
        for i,d in enumerate(dlist):
            sdata[d['id']]['name'] = 'HERV9_%s%s' % (band,someletters[i])
            if d['alias'] != '' and d['alias'] != '%s%s' % (band,someletters[i]):
                print >>sys.stderr, 'Mismatch: %s %s%s' % (d['alias'],band,someletters[i])
                d['alias'] += '*'                    

for mname,d in sdata.iteritems():
    if 'name' not in d:
        d['name'] = ''

columns =  ['id','category','chrom','band','alias','name']
print >>sys.stdout, '\t'.join(columns)
for mname in sorted(sdata.keys()):
    d = sdata[mname]
    print >>sys.stdout, '\t'.join([str(d[k]) for k in columns])

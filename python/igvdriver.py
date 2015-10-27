#! /usr/bin/env python

import os
from IGV import *
import re
import time

igv = IGV()
igv.new()
igv.genome('hg19')
igv.load(os.path.join(os.getcwd(),'../rmsk_LTR.hg19.gtf'))
igv.load(os.path.join(os.getcwd(),'concat.gtf'))
igv.load(os.path.join(os.getcwd(),'canonical.gtf'))
igv.load(os.path.join(os.getcwd(),'oneside.gtf'))
igv.load(os.path.join(os.getcwd(),'soloint.gtf'))
igv.load(os.path.join(os.getcwd(),'sololtr.gtf'))
igv.load(os.path.join(os.getcwd(),'unusual.gtf'))
igv.load(os.path.join(os.getcwd(),'subtables.gtf'))

# categories = ['oneside','soloint','unusual']
# for cat in categories:
#     if not os.path.exists(os.path.join(os.getcwd(),'%s_snapshots' % cat)):
#         os.mkdir(os.path.join(os.getcwd(),'%s_snapshots' % cat))
#     
#     igv.snapshotDirectory(os.path.join(os.getcwd(),'%s_snapshots' % cat))
#     lines = [l.strip('\n').split('\t') for l in open('transcript.gtf','rU')]
#     for l in lines:
#         if l[1] == 'merged':
#             attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
#             if attr['category'] == cat:
#                 print attr['name']
#                 locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
#                 igv.goto(locus)
#                 igv.expand()
#                 igv.snapshot(filename='%s.png' % attr['name'])

##########################################################################################

if not os.path.exists(os.path.join(os.getcwd(),'sub_snapshots')):
    os.mkdir(os.path.join(os.getcwd(),'sub_snapshots'))

igv.snapshotDirectory(os.path.join(os.getcwd(),'sub_snapshots'))

lines = [l.strip('\n').split('\t') for l in open('subtables.gtf','rU')]
for l in lines:
    locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
    attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
    print attr['locus']
    igv.goto(locus)
    igv.expand()
    igv.snapshot(filename='%s.png' % attr['locus'])




"""

### Oneside
if not os.path.exists(os.path.join(os.getcwd(),'oneside_snapshots')):
    os.mkdir(os.path.join(os.getcwd(),'oneside_snapshots'))

igv.snapshotDirectory(os.path.join(os.getcwd(),'oneside_snapshots'))

lines = [l.strip('\n').split('\t') for l in open('transcript.gtf','rU')]
for l in lines:
    if l[1] == 'merged':
        attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        if attr['category'] == 'oneside':
            print attr['name']
            locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
            igv.goto(locus)
            igv.expand()
            igv.snapshot(filename='%s.png' % attr['locus'])

### Internal
if not os.path.exists(os.path.join(os.getcwd(),'soloint_snapshots')):
    os.mkdir(os.path.join(os.getcwd(),'soloint_snapshots'))

igv.snapshotDirectory(os.path.join(os.getcwd(),'soloint_snapshots'))

lines = [l.strip('\n').split('\t') for l in open('transcript.gtf','rU')]
for l in lines:
    if l[1] == 'merged':
        attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        if attr['category'] == 'soloint':
            print attr['name']
            locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
            igv.goto(locus)
            igv.expand()
            igv.snapshot(filename='%s.png' % attr['locus'])


### Unusual
if not os.path.exists(os.path.join(os.getcwd(),'unusual_snapshots')):
    os.mkdir(os.path.join(os.getcwd(),'unusual_snapshots'))

igv.snapshotDirectory(os.path.join(os.getcwd(),'unusual_snapshots'))

lines = [l.strip('\n').split('\t') for l in open('transcript.gtf','rU')]
for l in lines:
    if l[1] == 'merged':
        attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
        if attr['category'] == 'unusual':
            print attr['name']
            locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
            igv.goto(locus)
            igv.expand()
            igv.snapshot(filename='%s.png' % attr['locus'])


"""
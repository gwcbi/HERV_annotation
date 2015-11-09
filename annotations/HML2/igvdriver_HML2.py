#! /usr/bin/env python

import os
import sys
from IGV import *
import re
import time

igv = IGV()
igv.new()
igv.genome('hg19')
igv.load(os.path.join(os.getcwd(),'../other_sources/rmsk_LTR.hg19.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/concat.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/prototype.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/oneside.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/soloint.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/sololtr.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/unusual.gtf'))
igv.load(os.path.join(os.getcwd(),'tmp/subtables.gtf'))

if not os.path.exists(os.path.join(os.getcwd(),'snapshots')):
    os.mkdir(os.path.join(os.getcwd(),'snapshots'))

categories = ['prototype', 'oneside', 'soloint', 'unusual']
lines = [l.strip('\n').split('\t') for l in open('transcript.gtf','rU')]
for cat in categories:
    if not os.path.exists(os.path.join(os.getcwd(),'snapshots/%s' % cat)):
        os.mkdir(os.path.join(os.getcwd(),'snapshots/%s' % cat))
    
    igv.snapshotDirectory(os.path.join(os.getcwd(),'snapshots/%s' % cat))
    for l in lines:
        if l[1] == 'merged':
            attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
            if attr['category'] == cat:
                print attr['name']
                locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
                igv.goto(locus)
                igv.expand()
                igv.snapshot(filename='%s.png' % attr['name'])

##########################################################################################

if not os.path.exists(os.path.join(os.getcwd(),'snapshots/literature')):
    os.mkdir(os.path.join(os.getcwd(),'snapshots/literature'))

igv.snapshotDirectory(os.path.join(os.getcwd(),'snapshots/literature'))

lines = [l.strip('\n').split('\t') for l in open('subtables.gtf','rU')]
for l in lines:
    locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
    attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
    print attr['locus']
    igv.goto(locus)
    igv.expand()
    igv.snapshot(filename='%s.png' % attr['locus'])

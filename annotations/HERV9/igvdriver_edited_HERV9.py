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
igv.load(os.path.join(os.getcwd(),'edited.hg19.gtf'))
igv.expand()

def get_coordinates_for_loci(loci,lines):
    for l in lines:
        if l[1] == 'merged':
            attr = dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8]))
            if attr['name'] in loci:
                locus = '%s:%d-%d' % (l[0],int(l[3])-5000, int(l[4])+5000)
                yield (attr['name'], locus)


merged = ['HERV9_2578','HERV9_0461','HERV9_0664','HERV9_1485','HERV9_1998','HERV9_2434','HERV9_3483','HERV9_3516']
joined = ['HERV9_0213','HERV9_1480','HERV9_1648','HERV9_2426','HERV9_2638','HERV9_3560','HERV9_4082','HERV9_4202','HERV9_4273']
lines = [l.strip('\n').split('\t') for l in open('edited.hg19.gtf','rU')]

for lname,lcoord in get_coordinates_for_loci(merged,lines):
    print lname
    igv.goto(lcoord)
    time.sleep(5)

for lname,lcoord in get_coordinates_for_loci(joined,lines):
    print lname
    igv.goto(lcoord)
    time.sleep(5)
#! /bin/bash

[[ ! -e initial_merge.hg19.gtf ]] && echo "ERROR: initial_merge.hg19.gtf not found" && exit 1

### Visualize current annotations ########################################################
mkdir -p tmp

cat ucsc/*.hg19.gtf | sortgtf > tmp/concat.gtf

cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'prototype' > tmp/prototype.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'oneside' > tmp/oneside.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'soloint' > tmp/soloint.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'sololtr' > tmp/sololtr.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'unusual' > tmp/unusual.gtf

python igvdriver_HERV9.py

echo -e '*\n!.gitignore' > snapshots/.gitignore
##########################################################################################

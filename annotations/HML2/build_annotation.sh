#! /bin/bash

### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc
echo -e '*\n!.gitignore' > ucsc/.gitignore

# Internal
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK-int";' > ucsc/HERVK-int.hg19.txt

# LTR
models="LTR5 LTR5A LTR5B LTR5_Hs"
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*
##########################################################################################

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done
##########################################################################################

### Determine merge distance #############################################################
mergedist=$(cat ucsc/HERVK-int.hg19.gtf | calculate_gaplength)
##########################################################################################

### Initial merge using bedtools cluster #################################################
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HML2 > initial_merge.hg19.gtf
##########################################################################################

### Setup GTF for visualization ##########################################################
mkdir -p tmp

cat ucsc/*.hg19.gtf | sortgtf > tmp/concat.gtf

cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'prototype' > tmp/prototype.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'oneside' > tmp/oneside.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'soloint' > tmp/soloint.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'sololtr' > tmp/sololtr.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'unusual' > tmp/unusual.gtf

cat ../other_sources/subramanianT*.hg19.gtf | sed 's/^/chr/' | sortgtf > tmp/subtables.gtf
##########################################################################################

### Take the snapshots ###################################################################
[[ $1 == 'SNAPSHOT' ]] && python igvdriver_HML2.py
##########################################################################################

### Manual merge/split ###################################################################
cat initial_merge.hg19.gtf | \
    manual_merge --names HML2_0706,HML2_0707 --category prototype | \
    manual_merge --names HML2_1025,HML2_1026 --category prototype | \
    manual_merge --names HML2_1067,HML2_1068 --category prototype | \
    manual_split --name HML2_1101 --split 3 --newname HML2_1177,HML2_1101 --category sololtr,prototype > edited.hg19.gtf
##########################################################################################

### Filter by covered length #############################################################
filter_covlen --threshold 500 < edited.hg19.gtf > filtered.hg19.gtf
##########################################################################################

### Assign names to loci #################################################################
# Map locus IDs to cytogenic bands
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b ../other_sources/cytoband.gtf | \
    perl -lne '/^chr([XY\d]+)\s.*name "(\S+)".*gene_id "([\S\.]+)"/;print "$2\t$1$3"' > tmp/cyto_name_map.txt

# Map locus IDs to aliases from Subramanian et al.
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b tmp/subtables.gtf | \
    perl -lne '/name "(\S+)".*locus "([\S\.]+)"/;print "$1\t$2"' > tmp/sub_name_map.txt

# Assign names to loci that are not solo LTR and are over 500 bp
# Also for loci that are already named
python names_HML2.py > tmp/name_table.txt
##########################################################################################

### Add locus tag to GTF #################################################################
add_locus_tag --mapping tmp/name_table.txt < filtered.hg19.gtf > final_combined.hg19.gtf
##########################################################################################

### Create final annotation files ########################################################
grep 'merged' final_combined.hg19.gtf > final_merged.hg19.gtf
grep -v 'merged' final_combined.hg19.gtf > final.hg19.gtf
##########################################################################################

gtf2table final_merged.hg19.gtf > final_table.hg19.txt
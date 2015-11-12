#! /bin/bash

### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH

### Download RepeatMasker tracks from UCSC ###############################################
mkdir -p ucsc

# Internal
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERV9-int";' > ucsc/HERV9-int.hg19.txt

# LTR
models="LTR12 LTR12B LTR12C LTR12D LTR12E LTR12F LTR12_"
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### Determine merge distance #############################################################
mergedist=$(cat ucsc/HERV9-int.hg19.gtf | calculate_gaplength)

### Initial merge using bedtools cluster #################################################
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HERV9 > initial_merge.hg19.gtf

### Setup GTF for visualization ##########################################################
mkdir -p tmp

cat ucsc/*.hg19.gtf | sortgtf > tmp/concat.gtf

cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'prototype' > tmp/prototype.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'oneside' > tmp/oneside.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'soloint' > tmp/soloint.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'sololtr' > tmp/sololtr.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'unusual' > tmp/unusual.gtf

### Take the snapshots ###################################################################
[[ $1 == 'SNAPSHOT' ]] && python igvdriver_HERV9.py

### Manual merge/split ###################################################################
cat initial_merge.hg19.gtf | \
        manual_merge --names HERV9_2578,HERV9_2579 --category prototype | \
        manual_merge --names HERV9_0461,HERV9_0462 --category prototype | \
        manual_merge --names HERV9_0664,HERV9_0665 --category prototype | \
        manual_merge --names HERV9_1485,HERV9_1486 --category prototype | \
        manual_merge --names HERV9_1998,HERV9_1999 --category prototype | \
        manual_merge --names HERV9_2434,HERV9_2435 --category prototype | \
        manual_merge --names HERV9_3483,HERV9_3484 --category prototype | \
        manual_merge --names HERV9_3516,HERV9_3517 --category prototype | \
        manual_split --name HERV9_0213 --split 5 --newname HERV9_0213,HERV9_4699 --category prototype,prototype | \
        manual_split --name HERV9_1480 --split 3 --newname HERV9_1480,HERV9_4700 --category prototype,soloint | \
        manual_split --name HERV9_1648 --split 1 --newname HERV9_4701,HERV9_1648 --category sololtr,prototype | \
        manual_split --name HERV9_2426 --split 4 --newname HERV9_2426,HERV9_4702 --category prototype,sololtr | \
        manual_split --name HERV9_2638 --split 1 --newname HERV9_4703,HERV9_2638 --category sololtr,prototype | \
        manual_split --name HERV9_3560 --split 2 --newname HERV9_3560,HERV9_4704 --category oneside,sololtr | \
        manual_split --name HERV9_4082 --split 2 --newname HERV9_4705,HERV9_4082 --category oneside,prototype | \
        manual_split --name HERV9_4202 --split 4 --newname HERV9_4202,HERV9_4706 --category oneside,sololtr | \
        manual_split --name HERV9_4273 --split 2 --newname HERV9_4273,HERV9_4707 --category oneside,sololtr > edited.hg19.gtf

### View edited files ####################################################################
[[ $1 == 'SNAPSHOT' ]] && python igvdriver_edited_HERV9.py

### Filter by covered length #############################################################
filter_covlen --threshold 500 < edited.hg19.gtf > filtered.hg19.gtf

### Assign names to loci #################################################################
# Map locus IDs to cytogenic bands
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b ../other_sources/cytoband.gtf | \
    perl -lne '/^chr([XY\d]+)\s.*name "(\S+)".*gene_id "([\S\.]+)"/;print "$2\t$1$3"' > tmp/cyto_name_map.txt

# Assign names to loci that are not solo LTR and are over 500 bp
python names_HERV9.py > tmp/name_table.txt

### Add locus tag to GTF #################################################################
add_locus_tag --mapping tmp/name_table.txt < filtered.hg19.gtf > HERV9_combined.hg19.gtf

### Create final annotation files ########################################################
grep 'merged' HERV9_combined.hg19.gtf > HERV9_merged.hg19.gtf
grep -v 'merged' HERV9_combined.hg19.gtf > HERV9.hg19.gtf
gtf2table HERV9_merged.hg19.gtf | add_subfamily --gtf  HERV9.hg19.gtf | HERV9.locus_table.txt

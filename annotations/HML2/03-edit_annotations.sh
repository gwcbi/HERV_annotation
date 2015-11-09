#! /bin/bash

### Check for initial_merge file #########################################################
[[ ! -e initial_merge.hg19.gtf ]] && echo "ERROR: initial_merge.hg19.gtf not found" && exit 1

### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH


### Determine annotations to merge/split  ################################################
### Merge:
# HML2_0706 + HML2_0707
#     Matches 11q12.3
# HML2_1025 + HML2_1026
#     HML2_1025 is a full length LTR5B(-) that has a full length LTR5B(+) insertion
#     HML2_1026 is 2315 bp away from HML2_1029, and in the same orientation.
#     HML2_1026 also matches the annotation for 19q13.41
# HML2_1067 + HML2_1068
#     Matches 20q11.22

### Split:
# HML2_1101
#     This is an LTR5B provirus with an LTR5Hs locus nearby.
#     Created merged annotation HML2_1177

### Confirmed:
# HML2_0323
#     HERVK-int inserted into 3' LTR, same as 4q32.1
# HML2_0411
#     internal LTR5B, same as 6p22.1
# HML2_0456
#     tandem provirus, same as 7p22.1a + 7p22.1b
# HML2_0793
#     HERVK-int in front of LTR5B
##########################################################################################

cat initial_merge.hg19.gtf | \
    manual_merge --names HML2_0706,HML2_0707 --category prototype | \
    manual_merge --names HML2_1025,HML2_1026 --category prototype | \
    manual_merge --names HML2_1067,HML2_1068 --category prototype | \
    manual_split --name HML2_1101 --split 3 --newname HML2_1177,HML2_1101 --category sololtr,prototype > edited.hg19.gtf

### Filter by covered length #############################################################
# Exclude locus if total number of bases aligned to model is less than threshold
filter_covlen --threshold 500 < edited.hg19.gtf > filtered.hg19.gtf
##########################################################################################


### Assign names to loci #################################################################
# Map locus IDs to aliases from Subramanian et al.
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b tmp/subtables.gtf | \
    perl -lne '/name "(\S+)".*locus "([\S\.]+)"/;print "$1\t$2"' > tmp/sub_name_map.txt

# Map locus IDs to cytogenic bands
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b ../other_sources/cytoband.gtf | \
    perl -lne '/^chr([XY\d]+)\s.*name "(\S+)".*gene_id "([\S\.]+)"/;print "$2\t$1$3"' > tmp/cyto_name_map.txt

# Assign names to loci that are not solo LTR and are over 500 bp
# Also for loci that are already named
python names_HML2.py > tmp/name_table.txt
# Mismatch: 3p12.3 3p12.3b
# Mismatch: 3q21.2 3q21.2b
# Mismatch: 4p16.1a 4p16.1b
# Mismatch: 4p16.1b 4p16.1d
# Mismatch: 6p22.1 6p22.1b
# Mismatch: 7p22.1b 7p22.1
# Mismatch: 8p23.1a 8p23.1c
# Mismatch: 8p23.1b 8p23.1e
# Mismatch: 8p23.1c 8p23.1f
# Mismatch: 8p23.1d 8p23.1g
# Mismatch: 11p15.4 11p15.4a
# Mismatch: 19p12c 19p12b
##########################################################################################

### Add locus tag to GTF #################################################################
add_locus_tag --mapping tmp/name_table.txt < filtered.hg19.gtf > final_combined.hg19.gtf
grep 'merged' final_combined.hg19.gtf > final_merged.hg19.gtf
grep -v 'merged' final_combined.hg19.gtf > final.hg19.gtf
##########################################################################################

gtf2table final_merged.hg19.gtf > final_table.hg19.txt
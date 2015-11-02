#! /bin/bash

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

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
##########################################################################################

### Set environment variables ############################################################
export PATH=/Users/bendall/Projects/HERV_annotation/HERV_annotation.git/bin:$PATH
export PYTHONPATH=/Users/bendall/Projects/HERV_annotation/HERV_annotation.git/python:$PYTHONPATH

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done
##########################################################################################

### Determine merge distance #############################################################
# The calculate_gaplength script is for examining the structure of Repeatmasker annotations
# Sometimes, a single hit found by Repeatmasker has long internal gaps; these are included
# in UCSC as two different annotation regions. I can identify these instances since the 
# annotations are spatially adjacent and have the same alignment score. (Although it is
# also possible that two separate hits are adjacent and have the same score).
#
# For the initial merging, the merge distance is chosen so that all internal annotations
# belonging to the same hit are considered together; thus the merge distance is equal to
# the maximum gap length for the internal models.

mergedist=$(cat ucsc/HERVK-int.hg19.gtf | calculate_gaplength)

# chr4:165920483-165920578(+) --- 308 --- chr4:165920886-165921385(+)
# chr6:28652000-28655158(+) --- 975 --- chr6:28656133-28659787(+)
# chr6:57624432-57626916(+) --- 24 --- chr6:57626940-57627179(+)
# chr17:7961226-7961800(+) --- 313 --- chr17:7962113-7965219(+)
# chr19:20389202-20393442(+) --- 297 --- chr19:20393739-20395211(+)
# chr19:53862347-53863266(+) --- 1006 --- chr19:53864272-53867053(+)
# chr20:32716337-32716882(+) --- 305 --- chr20:32717187-32717401(+)
# chr5:46002817-46008123(-) --- 707 --- chr5:46008830-46009038(-)
# chr6:42862411-42868111(-) --- 311 --- chr6:42868422-42869549(-)
# chr8:17765604-17769882(-) --- 291 --- chr8:17770173-17772457(-)
# chr9:139675748-139675934(-) --- 305 --- chr9:139676239-139676332(-)
# chr9:139678300-139679465(-) --- 1149 --- chr9:139680614-139682433(-)
# chr11:62143016-62143301(-) --- 296 --- chr11:62143597-62148586(-)
# chr19:28129492-28132457(-) --- 311 --- chr19:28132768-28137361(-)
# chr19:53248275-53251592(-) --- 82 --- chr19:53251674-53252083(-)
# min gap length:  24
# mean gap length: 445
# max gap length:  1149

##########################################################################################


### Initial merge using bedtools cluster #################################################
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HML2 > initial_merge.hg19.gtf
# canonical:     60
# oneside:     51
# soloint:     9
# sololtr:     1051
# unusual:     5
##########################################################################################


### Visualize current annotations ########################################################
mkdir -p tmp
cat ../other_sources/subramanianT*.hg19.gtf | sed 's/^/chr/' | sortgtf > tmp/subtables.gtf
cat ucsc/*.hg19.gtf | sortgtf > tmp/concat.gtf

cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'canonical' > tmp/canonical.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'oneside' > tmp/oneside.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'soloint' > tmp/soloint.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'sololtr' > tmp/sololtr.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'unusual' > tmp/unusual.gtf

python igvdriver_HML2.py

##########################################################################################

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

cp initial_merge.hg19.gtf edited.hg19.gtf
# Manually perform the edits described above
chmod a-w edited.hg19.gtf

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

echo -e '*\n!.gitignore' > tmp/.gitignore
echo -e '*\n!.gitignore' > snapshots/.gitignore
echo -e '*\n!.gitignore' > ucsc/.gitignore


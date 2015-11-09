Building HML2 Annotation
------------------------

## 0. Setup environment

Set your PATH and PYTHONPATH so that you can find the scripts in this repository

```bash
### Set environment variables
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
```

## 1. Download RepeatMasker Tracks

In RepBase, HML2 internal regions are identified as "ERVK-int" and the LTR is the LTR5 subfamily.
Members of the LTR5 subfamily include LTR5, LTR5A, LTR5B, and LTR5_Hs. In order to build
the annotation, first we need to download the RepeatMasker tracks from UCSC.

```bash
### Download RepeatMasker tracks from UCSC
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

# Output:
#      256 ucsc/HERVK-int.hg19.txt
#       21 ucsc/LTR5.hg19.txt
#      266 ucsc/LTR5A.hg19.txt
#      432 ucsc/LTR5B.hg19.txt
#      646 ucsc/LTR5_Hs.hg19.txt
#     1621 total
```

## 2. Reformat tracks as GTF

The next step is to format the RepeatMasker tracks as GTF. I have provided a script, called
`table2gtf` that puts the tracks into the desired format.

```bash
### Format UCSC tables as GTF
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

```

## 3. Determine merge distance

Since proviruses are often composed of several RepeatMasker hits, we need to merge annotations 
that are _nearby_ in order to identify the complete provirus. The problem is that there is
no clear definition of _nearby_. Instead of choosing a completely arbitrary cutoff, we use
information from the RepeatMasker tracks to inform our choice of merge distance. For the
initial merging, the merge distance is chosen so that all internal annotations belonging to
the same hit are considered together; thus the merge distance is equal to the maximum gap length for the internal models.

The calculate_gaplength script is for examining the structure of RepeatMasker annotations.
Sometimes, a single hit found by RepeatMasker has long internal gaps; these are included in 
UCSC as two different annotation regions. We identify these instances by looking for hits 
that are are spatially adjacent and have the same alignment score. (Since it is possible 
that two separate hits are adjacent and have the same score, we reject any hits that are 
more than 10kb apart).

```bash
mergedist=$(cat ucsc/HERVK-int.hg19.gtf | calculate_gaplength)

# Output:
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
# min gap length:    24
# mean gap length:   445
# median gap length: 308
# max gap length:    1149
```

## 4. Perform the initial merge

First, `bedtools` is used to cluster nearby annotations according to the merge distance. These
clusters are passed to the `mergeclusters` script, which categorizes each cluster depending
on the arrangement of features within the cluster. First, each annotation within the cluster
is labeled as either "LTR" or "Internal". For example, the configuration for a cluster might
initially look like this:

`LTR--LTR--Int--Int--Int--LTR`

Next we collapse this by combining adjacent annotations of the same type:

`LTR--Int--LTR`

Finally, the category for the cluster is assigned as follows:

  + prototype: Internal region flanked by LTRs on both sides: `LTR--Int--LTR`
  + oneside: Internal region flanked by LTRs on only one side: `LTR--Int` or `Int--LTR`
  + soloint: Solo internal (`Int`)
  + sololtr: Solo LTR (`LTR`)
  + unusual: Anything that doesn't fit into the above categories

The output is a GTF file with two types of annotations. The original annotations are preserved
so that downstream software can be aware of gaps in the locus, while "merged" annotations 
span all annotations for the cluster. and are used for grouping together the annotations.

```bash
### Initial merge using bedtools cluster
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HML2 > initial_merge.hg19.gtf

# Output:
# prototype:     60
# oneside:     51
# soloint:     9
# sololtr:     1051
# unusual:     5
```

## 5. Visualize initial merge results

Visual inspection of the merged annotations can quickly reveal whether your merging was
successful. In order to do this quickly, I split up the merged GTF file by locus type
(prototype, oneside, etc.). Next I load all these files into IGV and use the script 
`igvdriver_HML2.py` to control IGV. The script goes through the annotation and takes a 
snapshot of each locus. The snapshots can be easily paged through using normal image
viewing software.

```bash
mkdir -p tmp

cat ucsc/*.hg19.gtf | sortgtf > tmp/concat.gtf

cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'prototype' > tmp/prototype.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'oneside' > tmp/oneside.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'soloint' > tmp/soloint.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'sololtr' > tmp/sololtr.gtf
cat initial_merge.hg19.gtf | grep -v 'merged' | grep 'unusual' > tmp/unusual.gtf

cat ../other_sources/subramanianT*.hg19.gtf | sed 's/^/chr/' | sortgtf > tmp/subtables.gtf

python igvdriver_HML2.py
```

## 6. Identify annotations to merge/split

After the visual inspection, it was clear that there were a few loci that needed to be merged
and one locus that needed to be split. I used annotations from the literature to help in these
decisions.

#### Merge
  + HML2_0706 + HML2_0707: Matches 11q12.3
  + HML2_1025 + HML2_1026:
HML2_1025 is a full length LTR5B(-) that has a full length LTR5B(+) insertion.
HML2_1026 is 2315 bp away from HML2_1029, and in the same orientation.
HML2_1026 also matches the annotation for 19q13.41
  + HML2_1067 + HML2_1068: Matches 20q11.22

#### Split

  + HML2_1101
  This is an LTR5B provirus with an LTR5Hs locus nearby.
  Create a new merged annotation called HML2_1177

### Unusual
I also double-checked some of the unusual annotations to see what is going on:

  + HML2_0323: HERVK-int inserted into 3' LTR, same as 4q32.1
  + HML2_0411: internal LTR5B, same as 6p22.1
  + HML2_0456: tandem provirus, same as 7p22.1a + 7p22.1b
  + HML2_0793: HERVK-int in front of LTR5B

## 7. Manual merge/split

I wrote two scripts, `manual_merge` and `manual_split`, to handle editing the annotations.
These scripts correctly update the merge line(s) for the locus and rename all the annotations
belonging to the locus.

```bash
cat initial_merge.hg19.gtf | \
    manual_merge --names HML2_0706,HML2_0707 --category prototype | \
    manual_merge --names HML2_1025,HML2_1026 --category prototype | \
    manual_merge --names HML2_1067,HML2_1068 --category prototype | \
    manual_split --name HML2_1101 --split 3 --newname HML2_1177,HML2_1101 --category sololtr,prototype > edited.hg19.gtf
```

## 8. Filter by covered length

Exclude locus if total number of bases aligned to model is less than threshold

```bash
filter_covlen --threshold 500 < edited.hg19.gtf > filtered.hg19.gtf
```

## 9. Assign names to loci

Create a text file mapping the locus ID to the cytogenetic band.

```bash
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b ../other_sources/cytoband.gtf | \
    perl -lne '/^chr([XY\d]+)\s.*name "(\S+)".*gene_id "([\S\.]+)"/;print "$2\t$1$3"' > tmp/cyto_name_map.txt
```

Also create a text file mapping the locus ID to annotations from the literature

```bash
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b tmp/subtables.gtf | \
    perl -lne '/name "(\S+)".*locus "([\S\.]+)"/;print "$1\t$2"' > tmp/sub_name_map.txt

```

The script `names_HML2.py` creates names for each locus using the cytogenetic band. If multiple
loci are present in the same band, a letter (a,b,c...) is added to the name. We compare the
names we generate to the names given in the literature.

```bash
python names_HML2.py > tmp/name_table.txt

# Output:
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
```

## 10. Update GTF with locus name

The locus names generated in the previous step are incorporated into the GTF file. The name
is put into the "locus" field, which is used by telescope.

```bash
add_locus_tag --mapping tmp/name_table.txt < filtered.hg19.gtf > final_combined.hg19.gtf
```

## 11. Create final annotation files

Create a file containing only the merged lines, `final_merged.hg19.gtf`. This is useful
for overlap testing or simulation.

```bash
grep 'merged' final_combined.hg19.gtf > final_merged.hg19.gtf
```

Create a file containing only the annotation lines, `final.hg19.gtf`. This is used by telescope.

```bash
grep -v 'merged' final_combined.hg19.gtf > final.hg19.gtf
```

Create a table for each merged line.

```bash
gtf2table final_merged.hg19.gtf > final_table.hg19.txt
```

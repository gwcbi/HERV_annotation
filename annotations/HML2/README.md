Building HML2 Annotation
------------------------

## 0. Setup environment

Set your PATH and PYTHONPATH so that you can find the scripts in this repository

```bash
### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
```

## 1. Download RepeatMasker Tracks

In RepBase, HML2 internal regions are identified as "ERVK-int" and the LTR is the LTR5 subfamily.
Members of the LTR5 subfamily include LTR5, LTR5A, LTR5B, and LTR5_Hs. In order to build
the annotation, first we need to download the RepeatMasker tracks from UCSC:

```bash
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
# Output:
#      256 ucsc/HERVK-int.hg19.txt
#       21 ucsc/LTR5.hg19.txt
#      266 ucsc/LTR5A.hg19.txt
#      432 ucsc/LTR5B.hg19.txt
#      646 ucsc/LTR5_Hs.hg19.txt
#     1621 total
##########################################################################################
```

## 2. Reformat tracks as GTF

The next step is to format the RepeatMasker tracks as GTF. I have provided a script, called
`table2gtf` that puts the tracks into the desired format.

```bash
### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done
##########################################################################################
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
##########################################################################################
```

## 4. Perform the initial merge

```bash
### Initial merge using bedtools cluster #################################################
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HML2 > initial_merge.hg19.gtf
# Output:
# prototype:     60
# oneside:     51
# soloint:     9
# sololtr:     1051
# unusual:     5
##########################################################################################
```
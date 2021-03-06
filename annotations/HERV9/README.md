Building HERV9 Annotation
------------------------

## 0. Setup environment

Set your PATH and PYTHONPATH so that you can find the scripts in this repository

```bash
### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
```

## 1. Download RepeatMasker Tracks

In RepBase, HERV9 internal regions are identified as "HERV9-int" and the LTR is the LTR12 subfamily.
Members of the LTR12 subfamily include LTR12, LTR12B, LTR12C, LTR12D, LTR12E, LTR12F, and LTR12_. In order to build
the annotation, first we need to download the RepeatMasker tracks from UCSC.

```bash
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
```

```
# Output:
#     1036 ucsc/HERV9-int.hg19.txt
#      775 ucsc/LTR12.hg19.txt
#      212 ucsc/LTR12B.hg19.txt
#     2741 ucsc/LTR12C.hg19.txt
#      490 ucsc/LTR12D.hg19.txt
#      116 ucsc/LTR12E.hg19.txt
#      637 ucsc/LTR12F.hg19.txt
#      556 ucsc/LTR12_.hg19.txt
#     6563 total
```

## 2. Reformat tracks as GTF

The next step is to format the RepeatMasker tracks as GTF. I have provided a script, called
`table2gtf` that puts the tracks into the desired format.

```bash
### Format UCSC tables as GTF ############################################################
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
### Determine merge distance #############################################################
mergedist=$(cat ucsc/HERV9-int.hg19.gtf | calculate_gaplength)
```

```
# Output:
# chr1:199060761-199062438(+) --- 267 --- chr1:199062705-199065431(+)
# chr2:63549895-63552630(+) --- 61 --- chr2:63552691-63552960(+)
# chr2:117044761-117045749(+) --- 897 --- chr2:117046646-117046714(+)
# chr3:49583504-49583659(+) --- 310 --- chr3:49583969-49584012(+)
# chr3:190903792-190906265(+) --- 314 --- chr3:190906579-190909379(+)
# chr4:59329878-59330200(+) --- 310 --- chr4:59330510-59331114(+)
# chr5:54908982-54909013(+) --- 1682 --- chr5:54910695-54910840(+)
# chr5:114326132-114327597(+) --- 301 --- chr5:114327898-114334859(+)
# chr5:129972877-129973748(+) --- 320 --- chr5:129974068-129977515(+)
# chr6:32461907-32462120(+) --- 1362 --- chr6:32463482-32466216(+)
# chr6:32505885-32506239(+) --- 298 --- chr6:32506537-32506916(+)
# chr7:82270001-82271529(+) --- 306 --- chr7:82271835-82274645(+)
# chr8:29163035-29163920(+) --- 72 --- chr8:29163992-29164384(+) --- 69 --- chr8:29164453-29167891(+)
# chr8:49152709-49153005(+) --- 312 --- chr8:49153317-49153907(+)
# chr12:82290756-82291885(+) --- 23 --- chr12:82291908-82298973(+)
# chr12:107604180-107605511(+) --- 70 --- chr12:107605581-107606408(+)
# chr15:45302993-45304836(+) --- 277 --- chr15:45305113-45306908(+)
# chr17:26044791-26045158(+) --- 245 --- chr17:26045403-26047489(+)
# chr17:52364025-52366642(+) --- 282 --- chr17:52366924-52370973(+)
# chr21:18669696-18677759(+) --- 751 --- chr21:18678510-18678784(+)
# chrY:16270675-16277069(+) --- 294 --- chrY:16277363-16279269(+)
# chrY:17095616-17095731(+) --- 295 --- chrY:17096026-17100185(+)
# chrY:18385185-18385652(+) --- 77 --- chrY:18385729-18388838(+)
# chrY:25083780-25086220(+) --- 291 --- chrY:25086511-25088862(+)
# chrY:26718470-26719967(+) --- 302 --- chrY:26720269-26722620(+)
# chr1:39443981-39444575(-) --- 307 --- chr1:39444882-39444948(-)
# chr1:175073037-175074709(-) --- 69 --- chr1:175074778-175076240(-)
# chr2:9035761-9036086(-) --- 303 --- chr2:9036389-9036665(-)
# chr3:38480596-38480971(-) --- 296 --- chr3:38481267-38481691(-)
# chr3:49259676-49260843(-) --- 292 --- chr3:49261135-49262616(-) --- 297 --- chr3:49262913-49262949(-)
# chr3:151383373-151386792(-) --- 112 --- chr3:151386904-151388249(-)
# chr4:76965968-76969404(-) --- 69 --- chr4:76969473-76970900(-)
# chr4:173291176-173291558(-) --- 244 --- chr4:173291802-173292079(-)
# chr5:32554581-32555193(-) --- 57 --- chr5:32555250-32558999(-)
# chr6:64037776-64041197(-) --- 71 --- chr6:64041268-64042696(-)
# chr6:115553483-115560380(-) --- 69 --- chr6:115560449-115561882(-)
# chr8:52203208-52210388(-) --- 300 --- chr8:52210688-52211556(-)
# chr10:96524281-96532539(-) --- 298 --- chr10:96532837-96532970(-)
# chr15:45154269-45154554(-) --- 305 --- chr15:45154859-45156090(-)
# chr18:48644483-48644914(-) --- 94 --- chr18:48645008-48646437(-)
# chr19:9426762-9430376(-) --- 263 --- chr19:9430639-9431595(-)
# chrX:70720476-70720630(-) --- 300 --- chrX:70720930-70721748(-)
# chrY:21263614-21263637(-) --- 29 --- chrY:21263666-21266717(-)
# chrY:27239767-27242118(-) --- 302 --- chrY:27242420-27244860(-)
# min gap length:    23
# mean gap length:   301
# median gap length: 295
# max gap length:    1682
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
### Initial merge using bedtools cluster #################################################
cat ucsc/*.hg19.gtf | sortgtf | bedtools cluster -d $mergedist -i -| mergeclusters --prefix HERV9 > initial_merge.hg19.gtf
```

```
# Output:
# prototype:     265
# oneside:     133
# soloint:     33
# sololtr:     4260
# unusual:     7
```

## 5. Visualize initial merge results

Visual inspection of the merged annotations can quickly reveal whether your merging was
successful. In order to do this quickly, I split up the merged GTF file by locus type
(prototype, oneside, etc.). Next I load all these files into IGV and use the script 
`igvdriver_HERV9.py` to control IGV. The script goes through the annotation and takes a 
snapshot of each locus. The snapshots can be easily paged through using normal image
viewing software.

```bash
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
```

## 6. Identify annotations to merge/split

After the visual inspection:

#### Merge
  + **HERV9_2578 + HERV9_2579, prototype** Several insertions of HERV17-int and LTR12F
  + **HERV9_0461 + HERV9_0462, prototype** Close together in same direction
  + **HERV9_0664 + HERV9_0665, prototype** MER9a3 and HERVK9-int are inserted into HERV9(LTR12) provirus
  + **HERV9_1485 + HERV9_1486, prototype** LTR2B and Harlequin-int inserted into HERV9 provirus  
  + **HERV9_1998 + HERV9_1999, prototype** HERV17-int inserted into HERV9(LTR12F) provirus    
  + **HERV9_2434 + HERV9_2435, prototype** LTR17 and HERV17-int inserted into HERV9(LTR12B) provirus
  + **HERV9_3483 + HERV9_3484, prototype** Both are LTR12D in same direction
  + **HERV9_3516 + HERV9_3517, prototype** HERV30-int and MER11B are inserted into HERV9 (LTR12F) provirus
   
#### Split
  + **HERV9_0213** HERV9(LTR12F) and HERV9(LTR12\_v) proviruses are 500bp apart.
Annotations 1-5 assigned to HERV9\_0213 (prototype) and annotations 6-11 are assigned to new locus HERV9\_4699 (prototype).
  + **HERV9_1480** HERV9(LTR12F) provirus is near HERV9-int locus.
Annotations 1-3 assigned to HERV9\_1480 (prototype) and annotation 4 is assigned to new locus HERV9\_4700 (soloint)
  + **HERV9_1648** HERV9(LTR12) provirus is near solo LTR12C.
Annotation 1 assigned to new locus HERV9\_4701 (sololtr) and annotations 2-6 assigned to HERV9\_1648 (prototype).
  + **HERV9_2426** HERV9(LTR12\_) provirus is near solo LTR12B.
Annotations 1-4 assigned to HERV9\_2426 (prototype) and annotation 5 assigned to new locus HERV9\_4702 (sololtr).
  + **HERV9_2638** HERV9(LTR12/D) provirus is near solo LTR12C.
Annotation 1 assigned to new locus HERV9\_4703 (sololtr) and annotation 2-4 assigned to HERV9\_2638 (prototype).
  + **HERV9_3560** LTR12C and LTR12C in opposite directions.
Annotations 1-2 assigned to HERV9\_3560 (oneside) and annotation 3 assigned to new locus HERV9\_4704 (sololtr).
  + **HERV9_4082** HERV9(LTR12C) provirus is near one-sided locus in opposite direction.
Annotations 1-2 assigned to new locus HERV9\_4705 (oneside) and annotations 3-5 assigned to HERV9\_4082 (prototype).
  + **HERV9_4202** LTR12F + HERV9-int in opposite direction to LTR12C.
Annotations 1-4 assigned to HERV9\_4202 (oneside) and annotation 5 assigned to new locus HERV9\_4706 (sololtr).
  + **HERV9_4273** LTR12C and LTR12E in opposite directions.
Annotations 1-2 assigned to HERV9\_4273 (oneside) and annotation 3 assigned to new locus HERV9\_4707 (sololtr).

### Unusual
I also double-checked some of the unusual annotations to see what is going on:
  + **HERV9_0096** LTR12D is inserted into HERV9(LTR12D) provirus in opposite direction
  + **HERV9_1554** LTR12C is inserted into HERV9(LTR12B) provirus in same direction
  + **HERV9_1917** LTR12D is inserted into HERV9(LTR12) provirus in opposite direction

## 7. Manual merge/split

I wrote two scripts, `manual_merge` and `manual_split`, to handle editing the annotations.
These scripts correctly update the merge line(s) for the locus and rename all the annotations
belonging to the locus.

```bash
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
```

I wrote a script that instructs IGV to display all the edited loci, `igvdriver_edited_HERV9.py`.

```bash
### View edited files ####################################################################
[[ $1 == 'SNAPSHOT' ]] && python igvdriver_edited_HERV9.py
```

## 8. Filter by covered length

Exclude locus if total number of bases aligned to model is less than threshold

```bash
### Filter by covered length #############################################################
filter_covlen --threshold 500 < edited.hg19.gtf > filtered.hg19.gtf
```

## 9. Assign names to loci

Create a text file mapping the locus ID to the cytogenetic band. The script `names_HERV9.py`
creates names for each locus using the cytogenetic band. If multiple
loci are present in the same band, a letter (a,b,c...) is added to the name. We compare the
names we generate to the names given in the literature.

```bash
### Assign names to loci #################################################################
# Map locus IDs to cytogenic bands
grep 'merged' filtered.hg19.gtf | bedtools intersect -wo -a - -b ../other_sources/cytoband.gtf | \
    perl -lne '/^chr([XY\d]+)\s.*name "(\S+)".*gene_id "([\S\.]+)"/;print "$2\t$1$3"' > tmp/cyto_name_map.txt

# Assign names to loci that are not solo LTR and are over 500 bp
python names_HERV9.py > tmp/name_table.txt
```

## 10. Update GTF with locus name

The locus names generated in the previous step are incorporated into the GTF file. The name
is put into the "locus" field, which is used by telescope.

```bash
### Add locus tag to GTF #################################################################
add_locus_tag --mapping tmp/name_table.txt < filtered.hg19.gtf > HERV9_combined.hg19.gtf
```

## 11. Create final annotation files

Create two final annotation files from the combined file. `HERV9_merged.hg19.gtf` contains
only the merged lines, and `HERV9.hg19.gtf` contains only the annotation lines. The latter
is used by telescope. Also create a table containing one row for each merged line, `HERV9.locus_table.txt`.

```bash
### Create final annotation files ########################################################
grep 'merged' HERV9_combined.hg19.gtf > HERV9_merged.hg19.gtf
grep -v 'merged' HERV9_combined.hg19.gtf > HERV9.hg19.gtf
gtf2table HERV9_merged.hg19.gtf > HERV9.locus_table.txt
```

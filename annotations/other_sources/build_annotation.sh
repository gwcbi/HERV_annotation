#! /bin/bash

### Set environment variables ############################################################
export PATH=$(dirname $(dirname $PWD))/bin:$PATH
export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH

### Download the RepeatMasker table for all LTR ##########################################
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repClass = "LTR";' | table2gtf_short > rmsk_LTR.hg19.gtf
##########################################################################################

### Download the cytogenetic band table ##################################################
curl \
  --data 'db=hg19&hgta_track=cytoBand&hgta_table=cytoBand&hgta_regionType=genome&hgta_outputType=gff&boolshad.sendToGalaxy=0&boolshad.sendToGreat=0&boolshad.sendToGenomeSpace=0&hgta_outFileName=""&hgta_compressType=none&hgta_doTopSubmit="get output"' \
  https://genome.ucsc.edu/cgi-bin/hgTables | sortgtf | cytoband.gtf
##########################################################################################

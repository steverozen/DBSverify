#!/bin/bash

# This script runs in the linux shell on the ICGC Cancer Collaboratory.
# It fetches tumor and normal BAM files using the "SCORE" client, and
# then runs sambamba on each tumor-normal pair to create corresponding
# miniBAMs for ater download.
#
# Input:
#
# 1. A CSV file, each line of which specifies a BED file to use and
# the "object ID"s and file names of the corresponding tumor and normal
# BAM files. The name of the file is in variable MASTER_LINK below,
# currently "collaboratory_bams_2021_07_16.csv"; this file is in
# the steverozen/DBSverify github repo in the data-raw folder.
#
# 2.One BED file for each row in the MASTER_LINK file; the BED file
# specifies the positons to snip out of the tumour and normal BAM files to
# create a miniBAM with only the regions of interest. (The BED file
# is not fully specified, just the "alliquot id" is provided. The code
# below addes a suffix to get the file name.)
#
# Output: Two miniBAMs for each row in the MASTER_LINK file, one created from
# the tumor BAM file and one created from the normal BAM file; each miniBAM
# contains the genomic regions specified in the BED file.
#
# Note: this script does not download the miniBAMS to a "local" computer.

#Pre-req:
#0) latest sambamba installed in the instance
#1) score-client is installed and configured properly to run in the instance
#2) ~/bed_folder is created and all the relevant bed files loaded
#3) collaboratory_bams_2021_07_16.csv is in the ~/
#
# collaboratory_bams_2021_07_16.csv example content
#aliquot_id,icgc_donor_id,T_Specimen ID,T_Object ID,T_File Name,N_Specimen ID,N_Object ID,N_File Name
#0009b464-b376-4fbc-8a56-da538269a02f,DO46416,SP101724,2c9d90f1-5c46-5f3f-8a7f-771b36414cf2,bf5b874240ec430239013f4292d4151c.bam,SP101728,9ec1f43b-379c-58cd-85fe-3d09439058f1,7441677a3aed0864a23b17b84c0a28c9.bam
#003819bc-c415-4e76-887c-931d60ed39e7,DO36062,SP79365,cf1f4e65-f643-5d30-a2e0-b06892e357d1,f71ad7f9141d02dffb79accc8d0ea77a.bam,SP79719,4b971fd2-8240-5115-ad4c-f65db3e2c1c6,1499f95a62bc7ca67c0a6156e5189063.bam
#0040b1b6-b07a-4b6e-90ef-133523eaf412,DO45049,SP98853,94bf5398-2773-57f2-be44-55196a4a75ac,608fababf81bd2eebbaaf89e835d112b.bam,SP111847,40e95369-586a-5fc5-a789-a8b6d0989152,c87b70b95fa868a227d2ba1c11a80139.bam
#00508f2b-36bf-44fc-b66b-97e1f3e40bfa,DO48578,SP106808,2bd8c64e-c5b4-53c9-9079-c938b8c5406b,b81ee896228e57e6b9a977db485694b8.bam,SP106810,74ab8219-2c44-5784-be82-55b398fffa6a,24cd86d85a0ff9c4caddb5d072485a5f.bam

HOME_LOC=/home/centos
SOFTWARE_LOC=$HOME_LOC/software_folder
export JAVA_HOME=/usr/lib/jvm/jre-11
export PATH="$JAVA_HOME/bin:$PATH"
export STORAGE_PROFILE=collab

SCORE_CLIENT=$SOFTWARE_LOC/score-client-5.3.0/bin/score-client
SAMBAMBA=$SOFTWARE_LOC/sambamba-0.8.0-linux-amd64-static

BED_LOC=$HOME_LOC/bed_folder
OUTPUT_LOC=$HOME_LOC/bamSlice_folder
[ ! -d "$OUTPUT_LOC" ] && mkdir "$OUTPUT_LOC"
BAM_LOC=$HOME_LOC/bam_folder
[ ! -d "$BAM_LOC" ] && mkdir "$BAM_LOC"

MASTER_LINK="collaboratory_bams_2021_07_16.csv"

DO_ID_LIST=($(awk -F "," 'NR>1 {print $2}' $HOME_LOC/$MASTER_LINK))
BED_ID_LIST=($(awk -F "," 'NR>1 {print $1}' $HOME_LOC/$MASTER_LINK))

T_OBJECT_ID_LIST=($(awk -F "," 'NR>1 {print $4}' $HOME_LOC/$MASTER_LINK))
T_SP_ID_LIST=($(awk -F "," 'NR>1 {print $3}' $HOME_LOC/$MASTER_LINK))
T_BAM_FILE=($(awk -F "," 'NR>1 {print $5}' $HOME_LOC/$MASTER_LINK))
N_OBJECT_ID_LIST=($(awk -F "," 'NR>1 {print $7}' $HOME_LOC/$MASTER_LINK))
N_SP_ID_LIST=($(awk -F "," 'NR>1 {print $6}' $HOME_LOC/$MASTER_LINK))
N_BAM_FILE=($(awk -F "," 'NR>1 {print $8}' $HOME_LOC/$MASTER_LINK))

COUNT=0

for DO_ID in "${DO_ID_LIST[@]}"
do

#Skip hashsum validation due to time to calculate for a large bam file.
printf "Sample: $DO_ID:${T_SP_ID_LIST[$COUNT]}\nDownloading object id: ${T_OBJECT_ID_LIST[$COUNT]} corresponding to ${T_BAM_FILE[$COUNT]} using Collaboratory cloud ...\n"
$SCORE_CLIENT --profile collab download --validate false --object-id ${T_OBJECT_ID_LIST[$COUNT]} --output-dir $BAM_LOC

printf "Done.\nBam slicing using $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed ...\n"
$SAMBAMBA slice -o $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs.bam -L $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed $BAM_LOC/${T_BAM_FILE[$COUNT]}
if [ $? != 0 ]; then
echo "Bam slice of $BAM_LOC/${T_BAM_FILE[$COUNT]} using $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed failed for some reason  ....  exiting script."
exit 1
fi

printf "Done.\nSorting and indexing $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs.bam ..."
$SAMBAMBA sort -t 4 -o $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs_srt.bam $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs.bam
if [ $? != 0 ]; then
echo "Sambamba sort of $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs.bam failed for some reason  ....  exiting script."
exit 1
fi
rm $OUTPUT_LOC/${DO_ID}_${T_SP_ID_LIST[$COUNT]}_dbs.bam

#Skip hashsum validation due to time to calculate for a large bam file.
printf "Done.\nSample:$DO_ID:${N_SP_ID_LIST[$COUNT]}\nDownloading object id: ${N_OBJECT_ID_LIST[$COUNT]} corresponding to ${N_BAM_FILE[$COUNT]} using Collaboratory cloud ...\n"
$SCORE_CLIENT --profile collab download --validate false --object-id ${N_OBJECT_ID_LIST[$COUNT]} --output-dir $BAM_LOC

printf "Done.\nBam slicing using $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed from upstream side ...\n"
$SAMBAMBA slice -o $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs.bam -L $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed $BAM_LOC/${N_BAM_FILE[$COUNT]}
if [ $? != 0 ]; then
echo "Bam slice of $BAM_LOC/${N_BAM_FILE[$COUNT]} using $BED_LOC/${BED_ID_LIST[$COUNT]}_merged_PCAWG_DBS.bed failed for some reason  ....  exiting script."
exit 1
fi

printf "Done.\nSorting and indexing $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs.bam ..."
$SAMBAMBA sort -t 4 -o $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs_srt.bam $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs.bam
if [ $? != 0 ]; then
echo "Sambamba sort of $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs.bam failed for some reason  ....  exiting script."
exit 1
fi
rm $OUTPUT_LOC/${DO_ID}_${N_SP_ID_LIST[$COUNT]}_dbs.bam

printf "Done.\nRemoving downloaded bam files to conserve space ...\n"
rm $BAM_LOC/${N_BAM_FILE[$COUNT]} $BAM_LOC/${T_BAM_FILE[$COUNT]} $BAM_LOC/${N_BAM_FILE[$COUNT]}.bai $BAM_LOC/${T_BAM_FILE[$COUNT]}.bai
((COUNT++))
done
printf "Done.\nExiting script."
exit 0


#!/bin/bash
# a simple script to count reads with identical ids in 2 files
# usage: ./count_fastq.sh DIRNAME
# args: DIRNAME -- folder to search files in
# filenames must be listed in a file "ids.txt" one per line; paired 
# filenames must follow each other, e.g.
# pair1_1
# pair1_2
# pair2_1
# pair2_2
# etc.
if [ -z "$1" ]; then
        echo "Please specify input directory!"
        name=`basename "$0"`
        echo "Usage: ${name} DIRNAME"
        exit 101
fi
DIRNAME=${1%/} 
if [ ! -d $DIRNAME ]; then
        echo "path $DIRNAME does not exist!"
        exit 101
fi
names=(`cat ids.txt`)
OUTFILENAME="stats_${DIRNAME}.lst"
echo -e "id1\tid2\ttot_reads_1\tidentical" > $OUTFILENAME 
for (( idx = 0; idx < ${#names[@]}/2; idx++ )); do
        i=$((2*$idx))
        j=$((2*$idx + 1))
        grep ^@SRR "${DIRNAME}/${names[$i]}.fastq" | cut -f 1 -d " " | cut -d "." -f 2 > tmp1.tmp
        grep ^@SRR "${DIRNAME}/${names[$j]}.fastq" | cut -f 1 -d " " | cut -d "." -f 2 > tmp2.tmp
        cnt1=(`wc -l < tmp1.tmp`)
        cnt2=(`awk 'FNR==NR{array[$0]; next} ($0 in array) { count++ } END { print count } ' tmp1.tmp tmp2.tmp`)
        echo -e "${names[$i]}\t${names[$j]}\t${cnt1}\t${cnt2}" >> $OUTFILENAME
        rm tmp1.tmp tmp2.tmp
done

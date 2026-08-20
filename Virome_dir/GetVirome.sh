#!/bin/sh
# requirements in paths: sra-tools, seqyclean, sortmerna, BBmap, kaiju.

A=$(pwd)

cat $1 | while read line; do
    mkdir wd_dir
    cd wd_dir
    mkdir srtRNA_wd_tmp
    fasterq-dump $line 
    seqyclean -1 "${line}_1.fastq" -2 "${line}_2.fastq" -qual -o $line -c virome_databases/UniVec -minlen 50
    echo "seqyclean done"
    sortmerna -ref virome_databases/UniVec -reads "${line}_PE1.fastq" -reads "${line}_PE2.fastq" -a 10 -other ./$line -out2 -fastx -workdir srtRNA_wd_tmp
    echo "sortmerna done"
    repair.sh in1="${line}_fwd.fq" in2="${line}_rev.fq" out1="${line}_fwd_fixed.fastq" out2="${line}_rev_fixed.fastq"
    echo "repair done"
    kaiju -E 10-5 -v -z 10 -t virome_databases/nodes.dmp -f virome_databases/kaiju_db_viruses.fmi -i "${line}_fwd_fixed.fastq" -j "${line}_rev_fixed.fastq" -o ../$line.kaiju
    echo "kaiju done"
    cd $A
    rm -rf wd_dir
    echo "Processing of $line completed."

done


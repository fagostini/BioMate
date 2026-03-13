#!/usr/local/bin bash
set -e
set -o pipefail

differences=0
for sample in `find temp/Demultiplexing -name "*.fastq.gz" | grep -v Undetermined`; do
    original=`echo $sample | sed 's/temp\//temp\/20260310_LM43899_0385_A12GGASZR5\//'`
    zcat "$original" | grep -v @ | grep -e G -e A -e T -e C | sort > temp/a
    zcat "$sample" | grep -v @ | grep -e G -e A -e T -e C | sort > temp/b
    if [[ `diff temp/a temp/b | wc -l` -gt 0 ]]; then
        echo "$sample"
        differences=$((differences+1))
    fi
    rm temp/a temp/b
done

if [[ $differences -gt 0 ]]; then exit 1; fi

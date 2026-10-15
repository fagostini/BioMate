#!/usr/local/bin bash
set -e

differences=0
flowcell_id=${FLOWCELL_ID:-20260310_LM43899_0385_A12GGASZR5}
for sample in `find temp/Demultiplexed -name "*.fastq.gz" | grep -v -e Undetermined -e _I1_ -e _I2_`; do
    original="temp/$flowcell_id/Demultiplexing/${sample#temp/Demultiplexed/}"
    zcat "$original" | grep -v @ | grep -e G -e A -e T -e C | sort > temp/a
    zcat "$sample" | grep -v @ | grep -e G -e A -e T -e C | sort > temp/b
    if [[ `diff temp/a temp/b | wc -l` -gt 0 ]]; then
        echo "$sample"
        differences=$((differences+1))
    fi
    rm temp/a temp/b
done

echo ""

if [[ $differences -gt 0 ]]; then exit 1; fi

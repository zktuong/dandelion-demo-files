#!/bin/bash
set -eo pipefail

for GEX in `cut -f 1 -d "-" barcodes.txt | sort | uniq`
do
    echo $GEX
    #translate GEX to the actual TRAB ID
    SAMPLE=`grep $GEX meta.csv | cut -f 1 -d ,`
    echo $SAMPLE
    sshpass -f ~/.sshpass rsync -P kp9@farm5-login:/lustre/scratch117/cellgen/team205/sharedData/cs42/all-dandelion-paper-map/$SAMPLE/all_contig_annotations.csv .
    if [ ! -f conga.csv ]
    then
        head -n 1 all_contig_annotations.csv > conga.csv
    fi
    #prepend the GEX ID to the barcodes as that's how cs42 named these
    sed "s/^/$GEX-/" all_contig_annotations.csv | sed "s/,\([ACTG]\)/,$GEX-\1/" | grep -f barcodes.txt >> conga.csv
    rm all_contig_annotations.csv
done

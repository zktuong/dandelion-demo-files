#!/bin/bash
set -eo pipefail

for GEX in FCAImmP7292028 FCAImmP7292029 FCAImmP7292030 FCAImmP7292031 FCAImmP7292032 FCAImmP7292033 FCAImmP7292034 FCAImmP7292035 FCAImmP7528294 FCAImmP7528296 FCAImmP7555856 FCAImmP7555857 FCAImmP7555858 FCAImmP7555859 FCAImmP7555860 FCAImmP7555861 FCAImmP7555862 FCAImmP7579224 FCAImmP7579225 FCAImmP7579226 FCAImmP7579228 FCAImmP7579230 FCAImmP7579231 FCAImmP7579232 FCAImmP7803016 FCAImmP7803017 FCAImmP7803024 FCAImmP7803025 FCAImmP7803028 FCAImmP7803029 FCAImmP7803030 FCAImmP7803034 FCAImmP7803035 FCAImmP7851890 FCAImmP7851891 FCAImmP7851892 FCAImmP7851893 FCAImmP7851894 FCAImmP7851895 FCAImmP7964502 FCAImmP7964506 FCAImmP7964507 FCA_gut8015057 FCA_gut8015060 FCA_gut8015061 Human_colon_16S8159182 Human_colon_16S8159183 Human_colon_16S8159184 Human_colon_16S8159186 Human_colon_16S8159187 Human_colon_16S8159188 Human_colon_16S8159189 Human_colon_16S8159190
do
    echo $GEX
    #translate GEX to the actual TRAB ID
    SAMPLE=`grep $GEX meta.csv | cut -f 1 -d ,`
    echo $SAMPLE
    sshpass -f ~/.sshpass rsync -P kp9@farm5-login:/lustre/scratch117/cellgen/team205/sharedData/cs42/all-dandelion-paper-map/$SAMPLE/all_contig_annotations.csv .
    #prepend the GEX ID to the barcodes as that's how cs42 named these
    sed "s/^/$GEX-/" all_contig_annotations.csv | sed "s/,\([ACTG]\)/,$GEX-\1/" | grep -f barcodes.txt >> conga.csv
    rm all_contig_annotations.csv
done

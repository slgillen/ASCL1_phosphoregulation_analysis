#!/bin/bash

indexdir='indexes'

maindir=''
bwdir=$maindir/bigwig

conditions='d2_ESC d3_ESC d3_WT05 d3_SA01 d3_SA05'
for cond in ${conditions[@]}
do
    echo ${cond}
    wiggletools mean $bwdir/${cond}_*_RPGC.bw > $bwdir/${cond}_RPGC.bedgraph
    LC_COLLATE=C sort -k1,1 -k2,2n $bwdir/${cond}_RPGC.bedgraph > $bwdir/${cond}_RPGC_sorted.bedgraph
    bedGraphToBigWig $bwdir/${cond}_RPGC_sorted.bedgraph $indexdir/chromsizes.txt $bwdir/mean_${cond}_RPGC.bw
done
wait

#!/bin/bash

maindir='ATACseq'
datadir=$maindir/DiffBind
outdir=$maindir/motifanalysis_HOMER
mkdir $outdir


samplenames='' # add sample names for comparisons of interest

# run motif analysis for sites with increased accessibility
for s in $samplenames
do
    outdir=$maindir/${s}_inc_FDR005_sizegiven
    mkdir $outdir
    findMotifsGenome.pl $datadir/${s}_inc_FDR005.bed mm39 $outdir/ -size given -p 24 -mis 3 -S 20 
done
wait




#!/bin/bash

#mouse genome primary assembly: from https://www.gencodegenes.org/mouse/ version GRCm39 (mm39), release M29
#mouse gene annotation: from https://www.gencodegenes.org/mouse/ version GRCm39 (mm39), release M29

maindir=$1 
indexdir=$2 
samplefile=$3 
tempoutdir=$4 

indir=$maindir/raw_fastq 

#########  read in sample names from sample file  #########

echo "main directory is: $1"
echo "directory containing indexes is: $2"
echo "using sample sheet: $3"

while IFS=$'\r\n' read -r line 
do
	values=(${line//,/ })
	samplenames+=("${values[4]}") #condition and replicate info as sample name
	echo "${values[4]}"
done <<< $(tail +2 "$maindir/$samplefile")

echo "sample names are:" ${samplenames[@]}


######### initial Q.C check ########
fastqcdir=$maindir/initial_fastqc
mkdir $fastqcdir

for s in ${samplenames[@]}
do
  fastqc $indir/${s}_R1.fq.gz --outdir=$fastqcdir &
done
wait

for s in ${samplenames[@]}
do
  fastqc $indir/${s}_R2.fq.gz --outdir=$fastqcdir &
done
wait

echo "initial fastqc complete"


########  data filtering with fastp  ########
trimdir=$maindir/fastp
mkdir $trimdir

for s in ${samplenames[@]}
do	
  (fastp -i $indir/${s}_R1.fq.gz -I $indir/${s}_R2.fq.gz -w 4 -q 20 -l 20 -o $trimdir/${s}_fastp_R1.fq.gz -O $trimdir/${s}_fastp_R2.fq.gz) 2> $trimdir/fastp_stats_${s}.txt &
done
wait


########  Q.C check after trimming  ########
fastqcdir2=$maindir/fastqc_after_trimming
mkdir $fastqcdir2
for s in ${samplenames[@]}
do
  fastqc $trimdir/${s}_fastp_R1.fq.gz --outdir=$fastqcdir2 &
done
wait

for s in ${samplenames[@]}
do
  fastqc $trimdir/${s}_fastp_R2.fq.gz --outdir=$fastqcdir2 &
done
wait

echo "quality trimming and follow up QC complete"


# #########  align to genome with STAR aligner  ########
STARdir=$indexdir/STAR_indexes
aligndir=$tempoutdir/STARalignments
mkdir $aligndir

for s in ${samplenames[@]}
do
  STAR --readFilesCommand gunzip -c --readFilesIn $trimdir/${s}_fastp_R1.fq.gz $trimdir/${s}_fastp_R2.fq.gz --runThreadN 24 --genomeDir $STARdir --alignEndsType EndToEnd --outSAMtype BAM Unsorted --sjdbGTFfile $indexdir/gencode_vM29_annotation_ASCL1.gtf --outFileNamePrefix $aligndir/${s}
done
wait

echo "alignments complete"


########  sort the bam alignment files  ########
for s in ${samplenames[@]}
do	
  samtools sort -o $aligndir/${s}_sorted.bam $aligndir/${s}Aligned.out.bam &
done
wait

for s in ${samplenames[@]}
do	
  samtools flagstat $aligndir/${s}_sorted.bam > $aligndir/flagstat_${s}.txt &
done
wait


########  index the bam alignment files  ########
for s in ${samplenames[@]}
do
  samtools index $aligndir/${s}_sorted.bam $aligndir/${s}_sorted.bai & 
done
wait

for s in ${samplenames[@]}
do	
  rm $aligndir/${s}Aligned.out.bam &
done
wait

echo "alignments sorted and 
indexed"

# # ########  count the reads per gene  ########
aligndir=$tempoutdir/STARalignments
countsdir=$maindir/counts
mkdir $countsdir

for s in ${samplenames[@]}
do
  featureCounts -a $indexdir/gencode_vM29_annotation_ASCL1.gtf -s 2 -g gene_id --extraAttributes gene_name -B -C -p -o $countsdir/counts_${s}.txt $aligndir/${s}_sorted.bam &
done
wait


#######  format the counts table  ########
for s in ${samplenames[@]}
do
  cut -f1,7,8 $countsdir/counts_${s}.txt > $countsdir/counts_${s}_matrix.txt &
done
wait

#######  bring together all individual counts into one data table  ########
echo "Ensembl_id	gene_name	${samplenames[0]}" > $maindir/all_counts_WT_SA_ASCL1.txt

find "$countsdir/" -name "counts_${samplenames[0]}_matrix.txt" | xargs -I {} cat {} | sed -e '1,2d' >> $maindir/all_counts_WT_SA_ASCL1.txt

for s in ${samplenames[@]:1:${#samplenames[@]}}
do
	find "$countsdir/" -name "counts_${s}_matrix.txt" | xargs -I {} cat {} | sed -e '1,2d' | sed "1 i${s}" | cut -f3 > $maindir/all_counts_WT_SA_ASCL1_2.txt
	paste -d '\t' $maindir/all_counts_WT_SA_ASCL1.txt $maindir/all_counts_WT_SA_ASCL1_2.txt > $maindir/all_counts_WT_SA_ASCL1_3.txt
	mv $maindir/all_counts_WT_SA_ASCL1_3.txt $maindir/all_counts_WT_SA_ASCL1.txt
	rm $maindir/all_counts_WT_SA_ASCL1_2.txt
done

echo "counting complete"




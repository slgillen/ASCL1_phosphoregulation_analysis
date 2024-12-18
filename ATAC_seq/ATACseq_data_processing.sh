#!/bin/bash

#mouse genome primary assembly: from https://www.gencodegenes.org/mouse/ version GRCm39 (mm39), release M29
#mouse gene annotation: from https://www.gencodegenes.org/mouse/ version GRCm39 (mm39), release M29

maindir=$1
indexdir=$2 
samplefile=$3
indir=$maindir/raw_fastq 

#########  read in sample names from sample file  ########

echo "main directory is: $1"
echo "directory containing indexes is: $2"
echo "using sample sheet: $3"

while IFS=$'\r\n' read -r line 
do
	values=(${line//,/ })
	samplenames+=("${values[4]}") 
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


# ########  Q.C check after trimming  ########
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


########  alignment with bowtie2  ########
bt2indexdir=$indexdir/bowtie2_indexes
aligndir=$maindir/aligned
mkdir $aligndir


for s in ${samplenames[@]}
do
  echo ${s}
  (bowtie2 -p 40 -q -N 0 --dovetail --fr --no-mixed -x $bt2indexdir/GRCm39_primary_assembly_WT_ASCL1 -1 $trimdir/${s}_fastp_R1.fq.gz -2 $trimdir/${s}_fastp_R2.fq.gz --un-conc-gz $aligndir/unaligned_Mm_${s}.fq.gz -X 1000 --fr -S $aligndir/${s}.sam) 2> $aligndir/bowtie2_stats_${s}.txt
done
wait

########  sam to bam  ########
for s in ${samplenames[@]}
do
    samtools view -@ 20 -bS $aligndir/${s}.sam > $aligndir/${s}.bam
done
wait

for s in ${samplenames[@]}
do
    samtools flagstat --threads 20 $aligndir/${s}.bam > $aligndir/flagstat_${s}.txt &
done
wait

########  delete sam files  ########
for s in ${samplenames[@]}
do	
  rm $aligndir/${s}.sam &
done
wait

echo "to bam and flagstat complete"


########  filter bam  ########
filtdir=$maindir/filtering
mkdir $filtdir

for s in ${samplenames[@]}
do
    sambamba view -h -t 20 -f bam -F "[XS] == null and not unmapped and proper_pair" $aligndir/${s}.bam > $filtdir/${s}_filt.bam
done
wait

for s in ${samplenames[@]}
do
    samtools flagstat --threads 20 $filtdir/${s}_filt.bam > $filtdir/flagstat_${s}_filt.txt &
done
wait

echo "filter and flagstat complete"


########  sort the bam alignment files  ########

for s in ${samplenames[@]}
do	
  samtools sort -@ 20 -o $filtdir/${s}_filt_sorted.bam $filtdir/${s}_filt.bam 
done
wait


########  index the bam alignment files  ########
for s in ${samplenames[@]}
do
  samtools index $filtdir/${s}_filt_sorted.bam $filtdir/${s}_filt_sorted.bai & 
done
wait

echo "sort and index complete"


######## fragment size distribution ##########
QCdir=$maindir/ATAC_QC
mkdir $QCdir
for s in ${samplenames[@]}
do
    bamPEFragmentSize -p 20 -b $filtdir/${s}_filt_sorted.bam -hist $QCdir/PEfragsizes_${s}.png > $QCdir/fragsizeinfo_${s}.txt
done
wait


########  bam to bigwig  ########
bwdir=$maindir/bigwig
mkdir $bwdir

for s in ${samplenames[@]}
do
  bamCoverage -b $filtdir/${s}_filt_sorted.bam -o $bwdir/${s}_RPGC.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2728223909 --extendReads &
done
wait


########## peak calling ##########
regiondir=$maindir/regions
mkdir $regiondir

get the NFRs
for s in ${samplenames[@]}
do
  sambamba view -h -t 24 -f bam -F "((template_length > 0 and template_length < 135) or (template_length < 0 and template_length > -135))" $filtdir/${s}_filt.bam > $regiondir/NFR_${s}_sorted.bam
done
wait

for s in ${samplenames[@]}
do
  samtools index $regiondir/NFR_${s}_sorted.bam $regiondir/NFR_${s}_sorted.bai 
done
wait

peakdir=$maindir/peak_calling
mkdir $peakdir

for s in ${samplenames[@]}
do
  macs2 callpeak -t $regiondir/NFR_${s}_sorted.bam -g 2.72e9 -f BAMPE --keep-dup all --outdir $peakdir -n NFR_${s} -B 2> $peakdir/NFR_${s}_macs2.log &
done
wait


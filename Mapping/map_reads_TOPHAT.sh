##Default bash pipeline to map Ribo-seq, RNA-seq reads using Tophat2, some files could need specific pre-filtering (non-canonical adapters, UMIs, polyA...)
##Needed software: Tophat2, trim_galore, htseq-count

fastq=$1 #Files should end in '.fastq.gz'
gtf=$2
tophat_index=$3
contaminant_bowtie_index=$4
class=$5 #single,paired,ribo,rna29(RNA-seq reads trimmed to 29bp)
orient=$6 #yes,reverse, depending on the orientation of the reads for htseq-count

if [ $class == "single" ]
then
  bf1="$(basename -- $fastq)"
  echo "Filtering and trimming single RNA-seq "$bf1
  trim_galore --length 25 -q 30 --trim-n --gzip $fastq -o $bf1

  star="hg38_v98_index_rnaseq"
  ff1="${bf1//\.fastq.gz/_trimmed.fq.gz}"
  echo "Mapping "$ff1
  tophat --num-threads 6 --read-realign-edit-dist 0 -p 4 -z0 -M \ 
  -o $bf1 \
  -G $gtf \ 
  $tophat_index \
  $bf1/$ff1

  samtools index $bf1/accepted_hits.bam
  echo "Counting"
  htseq-count -f bam -r pos -s $orient $bf1/accepted_hits.bam $gtf > $bf1/htseq.exon.unique.counts
  htseq-count -f bam -r pos -s $orient -t CDS $bf1/accepted_hits.bam $gtf > $bf1/htseq.cds.unique.counts

elif [ $class == "paired" ]
then
  fastq1=$1
  fastq2=${fastq1/_1/_2}

  bf1="$(basename -- $fastq1)"
  echo "Filtering and trimming paired RNA-seq "$fastq1" "$fastq2
  trim_galore --paired --length 25 -q 30 --gzip $fastq1 $fastq2 -o $bf1

  ff1="${bf1//_1.fastq.gz/_1_val_1.fq.gz}"
  ff2="${bf1//_1.fastq.gz/_2_val_2.fq.gz}"
  echo "Mapping "$ff1" "$ff2
  tophat --num-threads 6 --prefilter-multihits --read-mismatches 2 --read-edit-dist 2 --read-realign-edit-dist 0 -z0 \ 
  -o $bf1 \
  -G $gtf \ 
  $tophat_index \
  $bf1/$ff1 \
  $bf1/$ff2

  samtools index $bf1/accepted_hits.bam
  echo "Counting"
  htseq-count -f bam -r pos -s $orient $bf1/accepted_hits.bam $gtf > $bf1/htseq.exon.unique.counts
  htseq-count -f bam -r pos -s $orient -t CDS $bf1/accepted_hits.bam $gtf > $bf1/htseq.cds.unique.counts

elif [ $class == "ribo" ]
then
  star="/fast/AG_Huebner/Jorge/mapping/hg38_v87_index_riboseq/"
  bf1="$(basename -- $fastq)"
  echo "Filtering and trimming single Ribo-seq "$bf1
  trim_galore --length 25 --trim-n --gzip $fastq -o $bf1

  star="hg38_v98_index_riboseq"
  ff1="${bf1//.fastq.gz/_trimmed.fq.gz}"
  echo "Filtering contaminants "$ff1
  bowtie2 --seedlen=25 -p 4 --un $bf1/$ff1\_filtered.fastq -x $contaminant_bowtie_index $bf1/$ff1>/dev/null

  echo "Mapping "$ff1
  tophat --num-threads 6 --read-realign-edit-dist 0 -p 4 -z0 -M \ 
  -o $bf1 \
  -G $gtf \ 
  $tophat_index \
  $bf1/$ff1

  samtools index $bf1/accepted_hits.bam
  echo "Counting"
  htseq-count -f bam -r pos -s $orient $bf1/accepted_hits.bam $gtf > $bf1/htseq.exon.unique.counts
  htseq-count -f bam -r pos -s $orient -t CDS $bf1/accepted_hits.bam $gtf > $bf1/htseq.cds.unique.counts

elif [ $class == "rna29" ]
then
  bf1="$(basename -- $fastq)"
  echo "Filtering and trimming single RNA-seq 29bp "$bf1
  trim_galore --length 29 --trim-n --gzip $fastq -o $bf1
  ff1="${bf1//.fastq.gz/_trimmed.fq.gz}"
  trim_galore --hardtrim5 29 --gzip $bf1/$ff1 -o $bf1

  star="hg38_v98_index_riboseq"
  ff2="${ff1//\.fq.gz/.29bp_5prime.fq.gz}"
  echo "Filtering contaminants "$ff2
  bowtie2 --seedlen=25 -p 4 --un $bf1/$ff2\_filtered.fastq -x $contaminant_bowtie_index $bf1/$ff2>/dev/null

  echo "Mapping "$ff2
  tophat --num-threads 6 --read-realign-edit-dist 0 -p 4 -z0 -M \ 
  -o $bf1 \
  -G $gtf \ 
  $tophat_index \
  $bf1/$ff2

  samtools index $bf1/accepted_hits.bam
  echo "Counting"
  htseq-count -f bam -r pos -s $orient $bf1/accepted_hits.bam $gtf > $bf1/htseq.exon.unique.counts
  htseq-count -f bam -r pos -s $orient -t CDS $bf1/accepted_hits.bam $gtf > $bf1/htseq.cds.unique.counts

fi


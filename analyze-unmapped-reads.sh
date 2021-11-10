#!/bin/bash

echo $0 $@

### USAGE

print_usage() {
  printf "

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ANALYZE UNMAPPED READS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 This script takes unmapped reads from aligned WGS data, or 
  paired fastq files generates scaffolds using the SPAdes 
  program, and blasts the resulting scaffolds against a 
  nucleotide database. 
 
 -------------------
 Required Programs
 -------------------

 This script uses the following programs:

  bedtools/2.27.1
  samtools/1.3.1
  blast/2.10.1+
  spades (https://github.com/ablab/spades)
  dustmasker (https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker)
  prinseq (http://prinseq.sourceforge.net)

 -------------------
 Flags
 -------------------

 -b (arg, required unless using -q)	BAM or CRAM file
 -o (arg, required)			the basename for the output file

 -a (arg)				number of blast alignments (default=5)
 -c (arg)				number of threads to allow (default=8)
 -d (arg)				dustmasker level (default=20)
 -f					only return unaligned fastq files
 -k					keep unmapped fastq files
 -l (arg)				path to the blast database (default=/n/shared_db/blastdb/202008)
 -m (arg)				min scaffold length to keep (default=500)
 -n (arg)				max percentage of Ns in order to keep a scaffold (default=10)
 -q (arg)				use pair1 fastq file as input
 -r					use RNA-Seq data (WARNING: Make sure your bam/cram file includes unmapped reads!)
 -s					skip blast step
 -t (arg)				a temporary directory to store intermediate files

 -h					print usage

 -------------------
 Example Usage
 -------------------

 analyze-unmapped-reads.sh -b sample.cram -o sample_ID

 -------------------
 Output
 -------------------

 outname-bam2fq-pair[12].fastq.gz	paired FASTQ files for the unmapped reads (requires -k flag)
 outname_spades_dmXX_ps.fasta		the list of scaffolds built from unmapped reads
 outname_blastn.txt			the blastn results for the scaffolds

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



"
}
### Get Arguments

cram=''
name=''
fq=''

nalign=5
threads=8
dm_level=20
fqonly=''
blastdb="/n/shared_db/blastdb/202008"
scamin=500
maxn=10
rna=''
skip=''
keep=''
tmp=.

spades=spades.py
dustmasker=dustmasker
prinseq=prinseq-lite.pl


while getopts 'b:o:a:c:d:fkl:m:n:t:q:rsh' flag; do
  case "${flag}" in
    b) cram="${OPTARG}" ;;
    o) name="${OPTARG}" ;;
    a) nalign="${OPTARG}" ;;
    c) threads="${OPTARG}" ;;
    d) dm_level="${OPTARG}" ;;
    f) fqonly="true" ;;
    k) keep="true" ;;
    l) blastdb="${OPTARG}" ;;
    m) scamin="${OPTARG}" ;;
    n) maxn="${OPTARG}" ;;
    t) tmp="${OPTARG}" ;;
    q) fq="${OPTARG}" ;;
    r) rna="true" ;;
    s) skip="true" ;;
    h) print_usage ; exit ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [ -z "$cram" ] && [ -z "$fq" ]; then
  print_usage
  echo "ERROR: -b and -q arguments missing, please enter a BAM/CRAM or fastq file"
  exit 1
elif [ -n "$cram" ] && [ -n "$fq" ]; then
  print_usage
  echo "ERROR: -b and -q flags are incompatible, please use one or the other"
  exit 1
elif [ -n "$fqonly" ] && [ -n "$fq" ]; then
  print_usage
  echo "ERROR: -f and -q flags are incompatible, please use one or the other"
  exit 1
elif [ -z "$name" ]; then
  print_usage
  echo "ERROR: -o argument missing, please enter an output name"
  exit 1
fi

if [ -n "$rna" ]; then
  echo "ATTN: Make sure your bam/cram includes unmapped reads (some RNA-Seq files exclude them)"
fi


### Setting variables

if [ -z "$skip" ]; then
  export BLASTDB=$blastdb
fi

outdir=${name}-out

tmp=${tmp}/${name}-temp
mkdir $tmp
mkdir $outdir

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Part 1: Get unmapped reads (samtools & bedtools)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Skip if using the -q flag
if [ -z "$fq" ]; then

  echo "Generating fastq files for unmapped reads"

  unmapped=${tmp}/${name}_unmapped.bam
  echo "Grabbing unmapped reads..."
  echo "This could take on the order of 1hr"
  date
  samtools view -@${threads} -f4 -F2048 -T $HG38 -bh $cram > $unmapped

  ## Get mates of unmapped reads (excluding unmapped reads themselves)

  mates_tmp=${tmp}/${name}_tmp_unmapped-mates.bam
  mates=${tmp}/${name}_unmapped-mates.bam

  echo "Grabbing mapped mates..."
  echo "This could take on the order of 1hr"
  date
  samtools view -@${threads} -f8 -F2048 -T $HG38 -bh $cram > $mates_tmp

  echo "Removing double counts (read & mate are both unmapped)..."
  date
  samtools view -@${threads} -F4 -T $HG38 -bh $mates_tmp > $mates

  rm $mates_tmp

  ## Merge unmapped and mates of unmapped reads

  merged=${tmp}/${name}_unmapped_mates.bam

  echo "Merging unmapped + mates..."
  date
  samtools merge -@${threads} $merged $unmapped $mates

  rm $unmapped $mates

  ## Sort merged file

  sorted=${tmp}/${name}_unmapped_mates.sorted.bam

  echo "Sorting..."
  date
  samtools sort -@16 -n $merged > $sorted

  rm $merged

  ## Convert bam to fastq

  r1=${outdir}/${name}-bam2fq-pair1.fastq
  r2=${outdir}/${name}-bam2fq-pair2.fastq

  echo "Converting to fastq..."
  date
  bamToFastq -i $sorted -fq $r1 -fq2 $r2

  bgzip $r1
  bgzip $r2

  r1=${r1}.gz
  r2=${r2}.gz

  if [ -n "$fqonly" ]; then
    rm -r $tmp
    echo "Unmapped reads can be found in $r1, $r2"
    exit 0
  fi

else

  r1=$fq
  r2=$(echo $r1 | sed 's/pair1/pair2/g' | sed 's/read1/read2/g')

  if [ -z "$(echo $r1 | grep 'pair1\|read1' | grep 'gz'$)" ]; then
    echo "ERROR: Please use zipped, paired fq files with the format \"pair1/2\" or \"read1/2\""
    exit 1
  fi

fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Part 2: Generate Scaffolds (spades)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

outdir_tmp=${tmp}/${name}_SPAdes

echo "Generating Scaffolds..."
date
$spades --only-assembler -m 22 -t ${threads} -1 $r1 -2 $r2 -o $outdir_tmp

scafa=${tmp}/scaffolds.fasta

cp ${outdir_tmp}/scaffolds.fasta $scafa

rm -r $outdir_tmp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Part 3: Clean Scaffolds (dustmasker & prinseq)
#
#  Removes short sequences and repeat sequences
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dm=${tmp}/${name}_dm${dm_level}.fa

echo "Masking repeats in scaffolds..."
date
$dustmasker -in $scafa -level $dm_level -outfmt fasta | sed 's/[cagtn]/N/g' | sed 's/leNNNh/length/g' | sed 's/Nov/cov/g' > $dm
rm $scafa

ps=${outdir}/${name}_spades_dm${dm_level}_ps

echo "Filtering out short scaffolds and repeat sequences..."
date
$prinseq -fasta $dm -out_format 1 -out_good $ps -min_len $scamin -ns_max_p $maxn -trim_tail_left 0 -trim_tail_right 0 -trim_ns_left 0 -trim_ns_right 0

fa=${ps}.fasta
rm $dm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Part 4: Blast Scaffold Sequences (blastn)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if [ -z "$rna" ] && [ -z "$skip" ] ; then
  blast_out=${outdir}/${name}_blastn.txt

  echo "Blasting scaffold sequences..."
  date
  blastn -num_threads ${threads} -outfmt "6 qseqid saccver evalue length sstart send ssciname scomname stitle sblastname sskingdom" -query $fa -out ${blast_out} -db nt -num_alignments $nalign
elif [ -z "$skip" ]; then
  blast_out=${outdir}/${name}_blastx.txt
  
  echo "Blasting scaffold sequences..."
  date
  blastx -num_threads ${threads} -outfmt "6 qseqid saccver evalue length sstart send ssciname scomname stitle sblastname sskingdom" -query $fa -out ${blast_out} -db refseq_protein -num_alignments $nalign
fi  

rm -r $tmp
if [ -z "$keep" ] && [ -z "$fq" ]; then
  rm $r1 $r2
fi

echo "Done!"
date








REF=$2
GTF=$3
RIBOBED=$4
BED=$5
PROJECT_DIR=$6
PBS_SCRIPT=$7

 
mkdir output_data 
mkdir bam
mkdir rseqc
mkdir fastqc_out
mkdir logfiles
mkdir featurecounts
cp $GTF.translation ./raw_data/gene_names.txt

sleep 5s

while read FASTQ1 FASTQ2 NAME Dx
do

echo $NAME

   qsub -v "R1=$PROJECT_DIR/fastq/$FASTQ1,R2=$PROJECT_DIR/fastq/$FASTQ2,OUTPUT=$PROJECT_DIR/bam/$NAME,REF=$REF,GTF=$GTF,RIBOBED=$RIBOBED,BED=$BED,POSTALN_OUT=$PROJECT_DIR/rseqc/$NAME,FEATURE_OUTPUT=$PROJECT_DIR/featurecounts/$NAME"  -N "$NAME" -o "$PROJECT_DIR/logfiles/rnaseq_pipeline_log_$NAME" $PBS_SCRIPT 
sleep 5s

done < $1



echo ---------------------------------------------------------------
PROJECT_DIR=/fs/scratch/osu0972/RNASeq5MS_5Healthy

echo Project directory:
echo $PROJECT_DIR 
echo

echo ---------------------------------------------------------------
PBS_SCRIPT=$PROJECT_DIR/hisat2_rnaseq_single.pbs
echo Project PBS script:
echo $PBS_SCRIPT
echo

echo ---------------------------------------------------------------
SAMPLELIST=$PROJECT_DIR/samples

echo Read sample info from:
echo $SAMPLELIST
echo "  Column 1: Fastq file list (R1)"
echo "  Column 2: Fastq file list (R2)"
echo "  Column 3: Sample base name"
echo

echo ---------------------------------------------------------------
REF=/fs/project/PCON0005/BISR/reference/HISAT2/Homo_sapiens/grch38/genome
     # e.g. REF=/fs/project/PCON0005/BISR/reference/GRCh38.primary/genome/star_index
echo Reference file:
echo $REF
echo

echo ---------------------------------------------------------------
GTF=/fs/project/PCON0005/BISR/reference/HISAT2/Homo_sapiens/grch38/Homo_sapiens.GRCh38.92.gtf
echo Gene Annotation File:
echo $GTF
echo

echo ---------------------------------------------------------------
RIBOBED=/fs/project/PCON0005/BISR/reference/HISAT2/Homo_sapiens/grch38/GRCh38_rRNA.bed
echo Ribosome File:
echo $RIBOBED
echo

echo ---------------------------------------------------------------
BED=/fs/project/PCON0005/BISR/reference/HISAT2/Homo_sapiens/grch38/hg38_GENCODE_v23.bed
echo Gene BED:
echo $BED
echo

echo $PROJECT_DIR > run.details
echo $PBS_SCRIPT >> run.details
echo $SAMPLELIST >> run.details
echo $REF >> run.details
echo $GTF >> run.details
echo $RIBOBED >> run.details
echo $BED >> run.details



./submit_pbs.sh $SAMPLELIST $REF $GTF $RIBOBED $BED $PROJECT_DIR $PBS_SCRIPT


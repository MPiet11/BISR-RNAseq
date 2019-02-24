#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -A PCON0005

cd $PBS_O_WORKDIR
module load R/3.5.0
export PATH=$PATH:/fs/project/PCON0005/BISR/programs/datamash-1.3
mkdir rseqc_reports
mkdir raw_data
mkdir output_data
 
echo "R1 R2 SampleID Base" | cat -   samples | awk '{print $3,$4}' > raw_data/sample.txt
sed -i 's/ /\t/g' raw_data/sample.txt

cd rseqc
#extracting insert size
files=(*.insert_size_metrics.txt) # store specific files in an array
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==7 {print $1,$2,$5,$6,$7,$8}'echo ${files[0]} | sed 's/^/sample_id\t/' > insert_size_metrics.txt
for i in *.insert_size_metrics.txt; do
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==8 {print ARGV[1],$1,$2,$5,$6,$7,$8}' $i >> insert_size_metrics.txt
done
sed -i 's/.insert_size_metrics.txt//g' insert_size_metrics.txt
mv insert_size_metrics.txt rseqc_insert_size_metrics.txt

#extracting mark duplicates
files=(*.marked_dup_metrics.txt) # store specific files in an array
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==7 {print $3,$9,$10}'echo ${files[0]} | sed 's/^/sample_id\t/' > duplication_metrics.txt
for i in *.marked_dup_metrics.txt; do
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==8 {print ARGV[1],$3,$9,$10}' $i >> duplication_metrics.txt
done
sed -i 's/.marked_dup_metrics.txt//g' duplication_metrics.txt
mv duplication_metrics.txt rseqc_duplication_metrics.txt

#extracting infer experiment
files=(*_infer) # store specific files in an array
awk 'BEGIN{FS=":"} FNR>=4 && FNR<=6 {print $1}'echo ${files[0]} | datamash transpose | sed 's/^/sample_id\t/' > Infer_metrics.txt
for i in *_infer; do
   awk 'BEGIN{FS=":"; print ARGV[1]} FNR>=4 && FNR<=6 {print $2}' $i | datamash transpose >> Infer_metrics.txt
done
sed -i 's/_infer/  /g' Infer_metrics.txt
sed -i 's/ /_/g' Infer_metrics.txt
sed -i 's/__//g' Infer_metrics.txt
sed -i 's/\t_/\t/g' Infer_metrics.txt
mv Infer_metrics.txt rseqc_Infer_metrics.txt

#extracting read distribution
files=(*_readdist) # store specific files in an array
sed -i 's/Total Reads/Total_Reads/' *_readdist
sed -i 's/Total Tags/Total_Tags/' *_readdist
sed -i 's/Total Assigned Tags/Total_Assigned_Tags/' *_readdist
awk 'FNR>=1 && FNR<=3 {print $1}'echo ${files[0]} | datamash transpose | sed 's/^/sample_id\t/' > read_distribution_P1.txt
for i in *_readdist; do
    awk 'BEGIN{print ARGV[1]} FNR>=1 && FNR<=3 {print $2}' $i | datamash transpose >> read_distribution_P1.txt
done
sed -i 's/_readdist//g' read_distribution_P1.txt
awk  -v OFS='\t' '{print $1,$2,$3,$4,$3-$4}'  read_distribution_P1.txt > read_distribution_P1_temp.txt
sed -i 's/\t0$/\tTotal_Unassigned_Tags/'  read_distribution_P1_temp.txt
mv  read_distribution_P1_temp.txt  read_distribution_P1.txt


awk 'FNR>=6 && FNR<=15 {print $1}'echo ${files[0]} | datamash transpose | sed 's/^/sample_id\t/' > read_distribution_P2.txt
for i in *_readdist; do
    awk 'BEGIN{print ARGV[1]} FNR>=6 && FNR<=15 {print $3}' $i | datamash transpose >> read_distribution_P2.txt
done
sed -i 's/_readdist//g' read_distribution_P2.txt

awk -v OFS='\t' 'FNR>=6 && FNR<=15 {print $1,$2}'echo ${files[0]} | datamash transpose > read_distribution_P3.txt
sed -i 's/^/Length\t/'  read_distribution_P3.txt
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$8,$11}' read_distribution_P3.txt > read_distribution_P3_temp.txt
mv read_distribution_P3_temp.txt read_distribution_P3.txt

paste  read_distribution_P1.txt read_distribution_P2.txt > read_distribution_P12.txt
awk -v OFS='\t' '{print $1,$2,$3,$4,$7,$8,$9,$10,$13,$16,$5}'  read_distribution_P12.txt > read_distribution_P12_temp.txt
mv read_distribution_P12_temp.txt read_distribution_P12.txt

mv read_distribution_P12.txt rseqc_read_distribution_P12.txt
mv read_distribution_P3.txt rseqc_read_distribution_P3.txt

mv rseqc_* ../rseqc_reports/

cd ../logfiles/
mv ../*.e* ./

grep Junc *.e* > temp_junc
sed -i 's/:Partial Novel Splicing Junctions:/\tPartial_Novel_Splicing_Junctions\t/' temp_junc
sed -i 's/:Novel Splicing Junctions:/\tNovel_Splicing_Junctions\t/' temp_junc
sed -i 's/:Known Splicing Junctions:/\tKnown_Splicing_Junctions\t/' temp_junc
sed -i 's/:Total splicing  Junctions:/\tTotal_Splicing_Junctions\t/' temp_junc
sed -i 's/\.e.......\t/\t/' temp_junc
awk   -v OFS='\t' '{print $1,$2,$3}' temp_junc | datamash -s crosstab 1,2 sum 3 > rseqc_junctions.txt

grep ribo rnaseq_* > temp_ribo
sed -i 's/rnaseq_pipeline_log_//' temp_ribo
sed -i 's/:/\t/' temp_ribo
sed -i 's/ (/\t/' temp_ribo
sed -i 's/):/\t/' temp_ribo
sed -i 's/ /_/g' temp_ribo
awk  -v OFS='\t' '{print $1,$3,$4}' temp_ribo | datamash -s crosstab 1,2 sum 3 > rseqc_ribo.txt
sed -i 's/Reads_consumed_by_input_gene_list/Ribosomal/' rseqc_ribo.txt
sed -i 's/Reads_not_consumed_by_input_gene_list/Not_ribosomal/' rseqc_ribo.txt
sed -i 's/qcfailed,_unmapped_reads/Unmapped/' rseqc_ribo.txt

mv rseqc_* ../rseqc_reports/

cd ../bam/
files=(*.hisat2.summary)
awk  -F':' 'FNR>=2 && FNR<=11 {print $1}'echo ${files[0]} | sed -e 's/\t\+//g' | datamash transpose | sed 's/^/sample_id\t/' > tmp_hisat2_metrics.txt
sed -i 's/ /_/g' tmp_hisat2_metrics.txt
for i in *.hisat2.summary; do
     awk -F':' 'BEGIN{print ARGV[1]} FNR>=2 && FNR<=11 {print $2}' $i | sed -e 's/\t\+//g' | datamash transpose  >> tmp_hisat2_metrics.txt
done
awk 'BEGIN {OFS=FS="\t"} {gsub(/\..*/,"",$1)}1' tmp_hisat2_metrics.txt | sed -e 's/([^(*)]*)//g' > hisat2_metrics.txt
sed -i 's/ //g' hisat2_metrics.txt
awk  -v OFS='\t'  '{print $1,$2,$6}' hisat2_metrics.txt > rseqc_hisat2_overall_alignment.txt
awk -v OFS='\t' '{print $1,$4,$5,$3}' hisat2_metrics.txt > rseqc_hisat2_uniq_multi_unmap.txt
 sed -i 's/Aligned_1_time\tAligned_>1_times\tAligned_0_time/UniquelyMapped\tMultiMapped\tUnMapped/'  rseqc_hisat2_uniq_multi_unmap.txt
mv hisat2_metrics.txt rseqc_hisat2_metrics.txt
mv rseqc_* ../rseqc_reports/

cd ../featurecounts/
grep -v Geneid *.featurecounts | grep -v '#' - | awk -v OFS='\t' '{print $1,$6,$7}' > temp_counts
sed -i 's/.featurecounts:/\t/' temp_counts
echo "gene_id length" > temp_length
awk '{print $2,$3}'  temp_counts | sort - | uniq > 1
cat 1 >> temp_length
sed -i 's/ /\t/' temp_length
awk -v OFS='\t' '{print $1,$2,$4}'  temp_counts | datamash -s crosstab 2,1 sum 3 > temp_counts_3

mv temp_length rseqc_featurecounts_length.txt
mv temp_counts_3 rseqc_featurecounts.txt

grep -v Status *.featurecounts.summary > temp_summary
sed -i 's/.featurecounts.summary:/\t/' temp_summary
datamash -s crosstab 2,1 sum 3 < temp_summary > rseqc_featurecounts_summary.txt

mv rseqc_* ../rseqc_reports/

cd ..
mkdir trash
mv bam/*.sam trash/
mv bam/*.bam trash/
mv trash/*sorted.bam bam/
mv bam/tmp_hisat2_metrics.txt trash/
#mv rseqc/*bed trash/
mv rseqc/*.r trash/
mv rseqc/*bam trash/
mv rseqc/*pdf trash/
mv rseqc/*xls trash/
mv featurecounts/temp* trash/

echo "if all files have been made successfully, you can delete everything in trash"

#mkdir raw_data
Rscript scripts/read_data_for_rshiny.R
#Rscript scripts/DGE_RNAseq_limma.R   
Rscript scripts/create_rshiny_input.R 

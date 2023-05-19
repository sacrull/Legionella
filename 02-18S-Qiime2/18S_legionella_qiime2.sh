#18S Legionella
#renaming
	#R1
ls *_R1_*.fastq | sort -n > num_order_R1_18S.txt
paste -d"\t" num_order_R1_18S.txt new_name_R1_18S.txt > name_R1_map_18S.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map_18S.txt
	#R2
ls *_R2_*.fastq | sort -n > num_order_R2_18S.txt
paste -d"\t" num_order_R2_18S.txt new_name_R2_18S.txt > name_R2_map_18S.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map_18S.txt	
	#moving
mv ./[mja]* ~/legionella/18S/data
gzip *.fastq

#manifest files header sample-id forward-absolute-filepath reverse-absolute-filepath
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manfist_18S.txt \
  --output-path ./paired-end-demux-18S.qza \
  --input-format PairedEndFastqManifestPhred33V2

cp paired-end-demux-18S.qza ~/legionella/18S


qiime demux summarize \
  --i-data ./paired-end-demux-18S.qza \
  --o-visualization ./paired-end-demux-sum-18S.qzv



#reverse compliments of primers used
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences paired-end-demux-18S.qza \
    --p-adapter-f AAAGGAATTGACGGARGGRCACC \
    --p-adapter-r CTACGACGGTATCTRATCRTCT \
    --p-quality-cutoff-3end 30 \
    --p-cores 4 \
    --verbose \
    --o-trimmed-sequences primer-trimmed-test3.qza
qiime demux summarize \
  --i-data ./primer-trimmed-test3.qza \
  --o-visualization ./primer-trimmed-test3-18S.qzv


#denoise data
qiime dada2 denoise-paired \
--i-demultiplexed-seqs primer-trimmed-test3.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 6 \
--o-representative-sequences representative-sequences-18S-test9.qza \
--o-table feature-table-18S-test9.qza \
--o-denoising-stats dada2-stats-18S-test9.qza

qiime metadata tabulate \
  --m-input-file dada2-stats-18S-test9.qza \
  --o-visualization dada2-stats-summ-18S-test9.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-18S-test9.qza \
  --o-visualization representative-sequences-summ-18S-test9.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-18S-test9.qza \
  --o-visualization representative-sequences-summ-18S-test9.qzv

biom convert -i feature-table.biom -o feature-table.txt --to-tsv


#silva

cp representative-sequences-18S.qza /home/suzanne/legionella/18S/silva
wget https://data.qiime2.org/2023.2/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2023.2/common/silva-138-99-tax.qza
	Michael S Robeson II, Devon R O’Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. RESCRIPt: Reproducible sequence taxonomy reference database management for the masses. bioRxiv 2020.10.05.326504; doi: https://doi.org/10.1101/2020.10.05.326504

qiime feature-classifier extract-reads \
    --i-sequences silva-138-99-seqs.qza \
    --p-f-primer AGAYGATYAGATACCGTCGTAG \
    --p-r-primer GGTGYCCYTCCGTCAATTCCTTT \
    --o-reads ref-seqs_SILVA-test1.qza

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads ref-seqs_SILVA-test1.qza \
 --i-reference-taxonomy silva-138-99-tax.qza \
 --o-classifier classifier_SILVA-test1.qza

qiime feature-classifier classify-sklearn \
  --i-classifier classifier_SILVA-test1.qza \
  --i-reads representative-sequences-18S-test9.qza \
  --p-n-jobs 6 \
  --o-classification taxonomy-test2.qza


qiime metadata tabulate \
  --m-input-file taxonomy-test2.qza \
  --o-visualization taxonomy-test2.qzv

qiime taxa filter-table \
  --i-table feature-table-18S-test9.qza \
  --i-taxonomy taxonomy-test2.qza \
  --p-mode contains \
  --p-include p__ \
  --p-exclude 'p__;,Unassigned,Mitochondria' \
  --o-filtered-table filtered-table1-18S-test1.qza


qiime feature-table filter-seqs \
  --i-data representative-sequences-18S-test9.qza \
  --i-table filtered-table1-18S-test1.qza \
  --o-filtered-data filtered-sequences1-18S-test1.qza


qiime taxa barplot \
  --i-table filtered-table1-18S-test1.qza \
  --i-taxonomy taxonomy-test2.qza \
  --m-metadata-file metadata_unsure_18S.txt \
  --o-visualization taxa-bar-plots-18S-test1.qzv


#pr2 v4.14.0



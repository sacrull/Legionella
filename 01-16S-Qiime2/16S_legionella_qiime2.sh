#16S Legionella
#renaming
	#R1
ls *_R1_*.fastq | sort -n > num_order_R1.txt
paste -d"\t" num_order_R1.txt new_name_R1.txt > name_R1_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map.txt
	#R2
ls *_R2_*.fastq | sort -n > num_order_R2.txt
paste -d"\t" num_order_R2.txt new_name_R2.txt > name_R2_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map.txt	
	#moving
mv ./[mja]* ~/legionella/16S/data
cp ./* ~/legionella/16S/
gzip *.fastq

#make qiime2 artifact
cd ~/legionella/16S/data

#manifest files header sample-id forward-absolute-path reverse-absolute-path
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manifest_16S_R.txt \
  --output-path ./paired-end-demux-16S-R.qza \
  --input-format PairedEndFastqManifestPhred33V2

cp paired-end-demux.qza ~/legionella/16S

qiime demux summarize \
  --i-data ./paired-end-demux-16S-R.qza \
  --o-visualization ./paired-end-demux-16S-R-vis.qza

#no cutadapt there is no primer

#denoise data
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux-16S-R.qza \
--p-trunc-len-f 240 \
--p-trunc-len-r 175 \
--p-n-threads 4 \
--o-representative-sequences representative-sequences-16S-R.qza \
--o-table feature-table-16S-R.qza \
--o-denoising-stats dada2-stats-16S-R.qza

qiime metadata tabulate \
  --m-input-file dada2-stats-16S-R.qza \
  --o-visualization dada2-stats-16S-summ.qzv
qiime feature-table summarize \
  --i-table feature-table-16S-R.qza \
  --m-sample-metadata-file metadata_unsure_16S-R.txt \
  --o-visualization feature-table-16S-R-summ.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-R.qza \
  --o-visualization representative-sequences-summ.qzv

biom convert -i feature-table.biom -o feature-table.txt --to-tsv

#taxanomy EzBioCloud
	#training classifier
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ezbiocloud_qiime_full.fasta \
--output-path ezbio-test1.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ezbiocloud_id_taxonomy.txt \
--output-path ref-taxonomy-test1.qza

qiime feature-classifier extract-reads \
	--i-sequences ezbio-test1.qza \
    --p-f-primer GTGCCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --o-reads ref-seqs_EZ_V3V4.qza

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads ref-seqs_EZ_V3V4.qza \
 --i-reference-taxonomy ref-taxonomy-test1.qza \
 --o-classifier classifier_EZ_V3V4-test1.qza
	actual classification
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_EZ_V3V4-test1.qza \
  --i-reads representative-sequences-R.qza \
  --p-n-jobs 4 \
  --o-classification taxonomy-test-16S-R.qza

qiime metadata tabulate \
  --m-input-file axonomy-test-16S-R.qza \
  --o-visualization taxonomy-test-R-vis.qzv


#taxonomy filter and rel abudance barchart
qiime taxa filter-table \
  --i-table ../data/feature-table-16S-R.qza \
  --i-taxonomy ../ezbio/taxonomy-test-16S-R.qza \
  --p-mode contains \
  --p-include "Bacteria;" \
  --o-filtered-table filtered-table-16S-R.qza


qiime feature-table filter-seqs \
  --i-data ../data/representative-sequences-16S-R.qza \
  --i-table filtered-table-16S-R.qza \
  --o-filtered-data filtered-sequences-16S-R.qza


qiime taxa barplot \
  --i-table ./filtered-table-16S-R.qza \
  --i-taxonomy ../ezbio/taxonomy-test-16S-R.qza \
  --m-metadata-file ../data/metadata_unsure_16S_R.txt \
  --o-visualization taxa-bar-plots-16S-filt-R.qzv 





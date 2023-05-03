# Qiime2 Pipeline for 16S reads
## 1. Renaming Fastq Files
My fastq files had names that I wanted to change to become more informative for the samples
### R1 Reads
Renaming R1 reads of fastq files
```
ls *_R1_*.fastq | sort -n > num_order_R1.txt
paste -d"\t" num_order_R1.txt new_name_R1.txt > name_R1_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map.txt
```
### R2 Reads
Renaming R2 reads of fastq files
```
ls *_R2_*.fastq | sort -n > num_order_R2.txt
paste -d"\t" num_order_R2.txt new_name_R2.txt > name_R2_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map.txt
```
## 2. Move files
Move files into a directory seperatre from the old fastq files
```
mv ./[mja]* ~/legionella/16S/data
cp ./* ~/legionella/16S/old
cd ~/legionella/16S/data
gzip *.fastq
```
## 3. Create Conda Enviroment
Use `qiime --help` to make sure Qiime2 is properly installed
```
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-linux-conda.yml
conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-linux-conda.yml
conda activate qiime2-2023.2
qiime --help
rm qiime2-2023.2-py38-linux-conda.yml
```
## 4. Make Qiime2 Artifact
Manifest should have sample-id forward-absolute-path reverse-absolute-path as headers
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manfist_16S.txt \
  --output-path ./paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
Move paired-end-demux.qza out of the current directory
```
cp paired-end-demux.qza ~/legionella/16S
```
Check the visualization to see where to trim
```
qiime demux summarize \
  --i-data ./paired-end-demux.qza \
  --o-visualization ./paired-end-demux-sum-viz.qzv
```
Export paired-end-demux-sum-viz.qzv and visualize it. I used cyberduck to export and used the qiime2 [web browser](https://view.qiime2.org/) for visualization
## 5. DADA2 pipeline
No cutadapt was used since there were no primers in the sequencing
### Make a directory
Make a directory for the denoising steps to occur in
```
mkdir denoise
cp paired-end-demux.qza ./denoise
cp metadata_16S.txt ./denoise
```
### Denoise data
```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f 240 \
--p-trunc-len-r 175 \
--p-n-threads 4 \
--o-representative-sequences representative-sequences-16S.qza \
--o-table feature-table-16S.qza \
--o-denoising-stats dada2-stats-16S.qza
```
### Visualize in order to check results
```
qiime metadata tabulate \
  --m-input-file dada2-stats-16S.qza \
  --o-visualization dada2-stats-16S-summ.qzv
qiime feature-table summarize \
  --i-table feature-table-16S.qza \
  --m-sample-metadata-file metadata_16S.txt \
  --o-visualization feature-table-16S-summ.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-16S.qza \
  --o-visualization representative-sequences-16S-summ.qzv
```
## 6. Assigning taxonomy
Ezbio cloud was used since my data was 16S
### Importing Ezbio Cloud
Request the database from their [website](https://www.ezbiocloud.net/resources/16s_download)
Upload the files into your directory
### Get all the files you will need into ./ezbio
Make sure you have the files from ezbio cloud
Get your files
```
cp ./representative-sequences-16S.qza ~/legionella/16S/ezbio/
cp
```
### Get files into Qiime2
```
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ezbiocloud_qiime_full.fasta \
--output-path ezbio-16S.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ezbiocloud_id_taxonomy.txt \
--output-path ref-taxonomy-16S.qza
```
### Get reads and train classifier
I used my primers for this part. You should use your own primers.
```
qiime feature-classifier extract-reads \
	--i-sequences ezbio-16S.qza \
    --p-f-primer GTGCCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --o-reads ref-seqs_EZ_V3V4.qza

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads ref-seqs_EZ_V3V4.qza \
 --i-reference-taxonomy ref-taxonomy-test1.qza \
 --o-classifier classifier_EZ_V3V4-test1.qza
	actual classification
```
### Assign taxonomy
```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_EZ_V3V4-test1.qza \
  --i-reads representative-sequences.qza \
  --p-n-jobs 4 \
  --o-classification taxonomy-16S.qza
```
Visualize in order to check taxonomy
```
qiime metadata tabulate \
  --m-input-file taxonomy-test1.qza \
  --o-visualization taxonomy-16S.qzv
```
## 7. Filter taxonomy
Only keep what was assigned to bacteria
```
qiime taxa filter-table \
  --i-table feature-table-0.qza \
  --i-taxonomy taxonomy-test1.qza \
  --p-mode contains \
  --p-include "Bacteria;" \
  --p-exclude d__Eukaryota;p__uncultured \
  --o-filtered-table filtered-table-16S.qza

qiime feature-table filter-seqs \
  --i-data representative-sequences-16S.qza \
  --i-table filtered-table-16S.qza \
  --o-filtered-data filtered-sequences-16S.qza
```
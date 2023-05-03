# Qiime2 Pipeline for 18S reads
## 1. Renaming Fastq Files
Move to your directory with your reads. I am renaming my files since I want more informative names. 
# R1 Reads
```
ls *_R1_*.fastq | sort -n > num_order_R1_18S.txt
paste -d"\t" num_order_R1_18S.txt new_name_R1_18S.txt > name_R1_map_18S.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map_18S.txt
```
# R2 Reads
```
ls *_R2_*.fastq | sort -n > num_order_R2_18S.txt
paste -d"\t" num_order_R2_18S.txt new_name_R2_18S.txt > name_R2_map_18S.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map_18S.txt
```
## 2. Moving files
Move renamed fastq files into new directory
```
mv ./[mja]* ~/legionella/18S/data
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
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manfist_18S.txt \
  --output-path ./paired-end-demux-18S.qza \
  --input-format PairedEndFastqManifestPhred33V2

cp paired-end-demux-18S.qza ~/legionella/18S

qiime demux summarize \
  --i-data ./paired-end-demux-18S.qza \
  --o-visualization ./paired-end-demux-sum-18S.qzv
```
## 5. Denoise
#Make denoise directory
```
cd ~/legionella/18S
mkdir denoise
cp paired-end-demux-18S.qza ./denoise
cd denoise
```
# Cutadapt
I had to remove the reverse compilment of primers since there was a lot of variety in my sequence read length
```
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences paired-end-demux-18S.qza \
    --p-adapter-f AAAGGAATTGACGGARGGRCACC \
    --p-adapter-r CTACGACGGTATCTRATCRTCT \
    --p-quality-cutoff-3end 30 \
    --p-cores 4 \
    --verbose \
    --o-trimmed-sequences primer-trimmed.qza
```
Visualize and check reads again
```
qiime demux summarize \
  --i-data ./primer-trimmed.qza \
  --o-visualization ./primer-trimmed-18S.qzv
```
# DADA2
No trimming was done because of the variable sequence length
```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs primer-trimmed.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 6 \
--o-representative-sequences representative-sequences-18S-.qza \
--o-table feature-table-18S.qza \
--o-denoising-stats dada2-stats-18S.qza
```
Visualize to make sure results look as expected and move representative-sequences-18S.qza into main 18S directory
```
qiime metadata tabulate \
  --m-input-file dada2-stats-18S.qza \
  --o-visualization dada2-stats-summ-18S.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-18S.qza \
  --o-visualization representative-sequences-summ-18S.qzv
qiime feature-table tabulate-seqs \
  --i-data representative-sequences-18S.qza \
  --o-visualization representative-sequences-summ-18S.qzv
cp representative-sequences-18S.qza ~/legionella/18S
```
## 6. Assign Taxonomy
Silva was used for 18S sequences

### Make silva directory
Get all files that will be needed into that directory
```
cd ~/legionella/18S
mkdir silva
cp representative-sequences-18S.qza /home/suzanne/legionella/18S/silva
wget https://data.qiime2.org/2023.2/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2023.2/common/silva-138-99-tax.qza
```
### Train classifier
```
qiime feature-classifier extract-reads \
    --i-sequences silva-138-99-seqs.qza \
    --p-f-primer AGAYGATYAGATACCGTCGTAG \
    --p-r-primer GGTGYCCYTCCGTCAATTCCTTT \
    --o-reads ref-seqs_SILVA-test1.qza

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads ref-seqs_SILVA.qza \
 --i-reference-taxonomy silva-138-99-tax.qza \
 --o-classifier classifier_SILVA-test1.qza
```
### Assign taxonomy
```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_SILVA.qza \
  --i-reads representative-sequences-18S.qza \
  --p-n-jobs 6 \
  --o-classification taxonomy-18S.qza
 ```

## 7. Filter taxonomy
Remove anything not assigned at phylum level or is assigned as mitochondrial
```
qiime taxa filter-table \
  --i-table feature-table-18S.qza \
  --i-taxonomy taxonomy-test2.qza \
  --p-mode contains \
  --p-include p__ \
  --p-exclude 'p__;,Unassigned,Mitochondria' \
  --o-filtered-table filtered-table1-18S.qza


qiime feature-table filter-seqs \
  --i-data representative-sequences-18S.qza \
  --i-table filtered-table1-18S.qza \
  --o-filtered-data filtered-sequences1-18S.qza
```






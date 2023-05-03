# Getting files out of Qiime2
```
conda activate qiime2-2023.2
```
## 16S
Get all files into same directory
### 1. Getting all files together
Get all files into same directory
```
mkdir qiime2R
cd ~/legionella/16S/denoise/
cp feature-table-16S.qza ~/legionella/R/qiime2R/
cp representative-sequences-16S.qza ~/legionella/R/qiime2R/
cd ~/legionella/16S/silva
cp taxonomy-16S.qza ~/legionella/R/qiime2R/
cd ~/legionella/R/qiime2R/
```
### 2. Feature table
#### 2.a Unzip feature table
```
unzip feature-table-16S
```
#### 2.b Convert from biom table
```
biom convert -i feature-table.biom -o feature-table1-16S.txt --to-tsv
```
#### 2.c Clean up
```
sed 's/# Constructed from biom file//' feature-table1-16S.txt | sed 's/#OTU ID/OTU/' | sed '1d' > feature-table-16S.txt 
```
#### 2.d Move to main R folder
```
cp feature-table-16S.txt ~/legionella/R
```
### 3. Taxnomy
#### 3.a Unzip taxonomy
```
unzip taxonomy-16S.qza
```
Go into the new directory and rename the file
```
cp taxonomy.tsv taxonomy-16S.tsv
```
Download file to computer. Get rid of the confidence value column before running the next steps.
#### 3.b Edit taxonomy
Makes taxonomy table have no empty columns and no spaces. Make sure the script and taxonomy file are in the same directory
```
python3 fix_taxonomy_16S.py taxonomy-16S.tsv > taxonomy_16S.txt

```
Manually look through taxonomy to make sure it looks normal.
#### 3.c Clean up taxonomy
Upload taxonomy back to machine into the `~/legionella/R`
Get rid of __ and taxa level indicators
```
sed 's/;/\t/g' taxonomy_16S.txt | sed '1d' > temp
mv temp taxonomy_18S.txt
```
### 4. Representative Sequences
#### 4.a Unzip sequences
```
unzip representative-sequences-16S.qza
```
#### 4.b Rename file
Go into the directory, rename the file, and move it.
```
cp dna-sequences.fasta ~/legionella/R/16S_rep_seqs.fasta
```

## 18S
Get all files into same directory
### 1. Getting all files together
Get all files into same directory
```
mkdir qiime2R
cd ~/legionella/18S/denoise/
cp feature-table-18S.qza ~/legionella/R/qiime2R/
cp representative-sequences-18S.qza ~/legionella/R/qiime2R/
cd ~/legionella/18S/silva
cp taxonomy-18S.qza ~/legionella/R/qiime2R/
cd ~/legionella/R/qiime2R/
```
### 2. Feature table
#### 2.a Unzip feature table
```
unzip feature-table-18S
```
#### 2.b Convert from biom table
```
biom convert -i feature-table.biom -o feature-table1-18S.txt --to-tsv
```
#### 2.c Clean up
```
sed 's/# Constructed from biom file//' feature-table1-18S.txt | sed 's/#OTU ID/OTU/' | sed '1d' > feature-table-18S.txt 
```
#### 2.d Move to main R folder
```
cp feature-table-18S.txt ~/legionella/R
```
### 3. Taxnomy
#### 3.a Unzip taxonomy
```
unzip taxonomy-18S.qza
```
Go into the new directory and rename the file
```
cp taxonomy.tsv taxonomy-18S.tsv
```
Download file to computer. Get rid of the confidence value column before running the next steps.
#### 3.b Edit taxonomy
Makes taxonomy table have no empty columns and no spaces. Make sure the script and taxonomy file are in the same directory
```
python3 fix_taxonomy_18S.py taxonomy-18S.tsv > taxonomy_18S.txt

```
Manually look through taxonomy to make sure it looks normal.
#### 3.c Clean up taxonomy
Upload taxonomy back to machine into the `~/legionella/R`
Get rid of __ and taxa level indicators
```
sed 's/d__//g' taxonomy_18S.txt | sed 's/;[pcfogs]__/\t/g' | sed '1d' > temp
mv temp taxonomy_18S.txt
```
### 4. Representative Sequences
#### 4.a Unzip sequences
```
unzip representative-sequences-18S.qza
```
#### 4.b Rename file
Go into the directory, rename the file, and move it.
```
cp dna-sequences.fasta ~/legionella/R/18S_rep_seqs_silva.fasta
```












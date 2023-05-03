Exiting Qiime2
#16S go under ezbio to find these
cp taxonomy-test-16S-R.qza ~/legionella/R/qiime2R/taxanomy-test1-16S.qza
unzip taxanomy-test1-16S.qza
#change using python
cp metadata_16S.txt ~/legionella/R/qiime2R
cp representative-sequences.qza ~/legionella/R/qiime2R/ref-seqs-16S.qza
cp feature-table-16S-R.qza ~/legionella/R/qiime2R/feature-table-16S.qza
unzip taxanomy-test1-16S.qza
biom convert -i feature-table.biom -o feature-table1-16S.txt --to-tsv
sed 's/# Constructed from biom file//' feature-table1-16S.txt | sed 's/#OTU ID/OTU/' | sed '1d' > feature-table-16S.txt

#18S check under silva to find
cp feature-table-18S-R.qza ~/legionella/R/qiime2R/
cp metadata_18S.txt ~/legionella/R/qiime2R
cp representative-sequences-18S-R.qza ~/legionella/R/qiime2R/
cp taxonomy-test-18S-R.qza ~/legionella/R/qiime2R/
unzip feature-table-18S-R.qza
biom convert -i feature-table.biom -o feature-table1-18S.txt --to-tsv
#18S sequence table
sed 's/# Constructed from biom file//' feature-table1-18S.txt | sed 's/#OTU ID/OTU/' | sed '1d' > feature-table-18S.txt 
unzip taxonomy-test-18S-R.qza
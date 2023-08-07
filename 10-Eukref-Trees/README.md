# Eukref Tree and Heatmap
Example shown is for Korotnevella tree and heatmap. All trees and heatmap follow the same method.
## 1. Install
Create conda envioment for tree
```bash
cd
conda create -n eukref
cd bin
wget https://drive5.com/downloads/usearch8.0.1623_i86linux32.gz
gzip -d usearch8.0.1623_i86linux32.gz
chmod +x usearch8.0.1623_i86linux32
mv usearch8.0.1623_i86linux32 usearch
wget http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz
tar -xvzf trimal.v1.2rev59.tar.gz
cd trimAl/source
make
mv trimal ../..
conda activae eukref
conda install -c conda-forge biopython
conda install -c bioconda seqtk
conda install -c "bioconda/label/cf201901" raxml
conda install -c "bioconda/label/cf201901" vsearch
conda install -c "bioconda/label/cf201901" mafft
conda install -c bioconda sina
#add to profile
export PATH="$PATH:/home/suzanne/bin"
```
Create conda enviroment for heatmap
```bash
#heatmap
conda deactivate
conda create -n phyloseq_tree
conda activate phyloseq_tree
conda install r-devtools
conda install -c conda-forge r-ggpubr
```
Open up R
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("reshape2") 
install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("ape")
install_github("microbiome/microbiome")
```
## 2.Backbone Tree Steps
### 2a. Make directory and get databases
Should have 3 databases total: pr2, silva, and NCBI
```bash
conda activate eukref
#mkdir
mkdir Korotnevella
#get databases
wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz
gzip -d *gz
#clean up silva databse
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' SILVA_138_SSURef_tax_silva.fasta > temp
mv temp SILVA_138_SSURef_tax_silva.fasta
#remove ward wrap for silva database
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' SILVA_138_SSURef_tax_silva.fasta > temp
mv temp SILVA_138_SSURef_tax_silva.fasta
#get korotnevella 
grep "Korotnevella" pr2_version_5.0.0_SSU_taxo_long.fasta -A 1 | sed 's/--//' > pr2_koro.fa
grep "Korotnevella" SILVA_138_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_koro.fa
ncbi query: Korotnevella 18S rRNA NOT whole genome shotgun NOT mrna 
sed '/^[^>]/s/U/T/g' silva_koro.fa > temp
mv temp silva_koro.fa
#remove word wrap for ncbi database and reformat header
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' sequence_koro.fasta > temp
mv temp ncbi_koro.fa
#get all refrences together
cat ncbi_koro.fa silva_koro.fa pr2_koro.fa > all_koro.fa
```
### 2b. Remove unwanted genomes
If you have created a tree and want to remove any outliers you can use this. Add unwanted genomes to remove.txt 
```bash
mkdir remove
nano remove.txt
cd remove
sed 's/_.*//g' ../all_koro.fa > temp
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' temp
for i in `cat ../remove.txt`; do rm ${i}.fasta; done
cat *.fasta > filt_koro.fa
cd ..
```
### 2c. Sort and Cluster Genomes using Vsearch
You can either use the removed genomes here or keep going with no genomes filtered out.
```bash
#sort by and cluster by 99%
vsearch --sortbylength ./all_koro.fa --output all_koro.sort #when going back change ./all_koro.fa to ../remove/filt_acan.fa
vsearch --cluster_fast all_koro.sort --centroids all_koro.clust --id 0.99
```
### 2d. Align using MAFFT
```bash
#clean up headers
sed 's/|.*//' all_koro.clust > temp
mv temp all_koro.clust
#using mafft
mafft --auto all_koro.clust > all_koro.align.fa
```
### 2e. Trim using Trimal
```bash
#clean up for trimal
sed '/^[^>]/s/n/-/g' all_koro.align.fa > temp
mv temp all_koro.align.fa
#use trimal
~/bin/trimal -in all_koro.align.fa -out all_koro.trim.fa  -gt 0.3 -st 0.001
```
### 2f. Backbone tree using RAxML
If you do not like how reference tree looks you can either add the genomes you want to remove to remove.txt or trim it using the archaeopteryx tree viewer.
```bash
#make backbone tree
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 7 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n koro_ref.tre -s all_koro.trim.fa
#organize
mkdir reference_trees
grep ">" all_koro.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.koro.txt
mv annotations.koro.txt reference_trees/
mv RAxML_bipartitions.koro_ref.tre reference_trees
cd ..
```
## 3. Make tree with your ASVs
### 3a. Get your ASVs
Use your taxonomy file to get their IDs and then pull their full fasta file
```bash
awk '{print $1, $7}' ../18S-tax-renamed.txt | grep Korotnevella| awk '{print $1}' > koro.ids
#pull their fasta file
cat koro.ids | while read line; do grep -w $line ../ref-seq-18S.fasta -A 1; done > koro.fa
```
### 3b. Align your ASVs with referene genomes
```bash
#allign sequences
mkdir placement_tree
cd placement_tree

sina -i ../all_koro.align.fa --prealigned -o all_koro.arb
sina -i ../koro.fa -r all_koro.arb -o query.align.koro.fa --fs-msc 0.01 --fs-full-len=100
```
### 3c. Concat reference sequences with your ASVs
```bash
#concat refrence sequences
cat query.align.koro.fa ../all_koro.align.fa > queryPlus.align.koro.fa
```
### 3d. Create tree
```bash
#create tree
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n koro.all.tre -s queryPlus.align.koro.fa -t ../reference_trees/RAxML_bipartitions.koro_ref.tre -T 7
#clean up tree
sed 's/QUERY___//g' RAxML_labelledTree.koro.all.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.koro.epa.tre
```
### 3e. Make constraint tree
Constraint tree will be used for alignment for heatmap. Visualize the tree using Figtree
```bash
#get constraint tree
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n koro.cons.tre -s queryPlus.align.koro.fa -g ../reference_trees/RAxML_bipartitions.koro_ref.tre -T 7 -x 25734 -p 25793
#export for midpoint rooting
RAxML_bipartitions.koro.cons.tre
```
### 3f. Annotations file
```bash
#make complete annotation file
awk '{print $1,$NF}' ../koro.ids | sed 's/ /\t/g' | sed 's/ASV/ASV /2' > my_koro_annots.txt
cat ../reference_trees/annotations.koro.txt my_koro_annots.txt > all_koro_annotations.txt
```
## 4. Heatmap
### 4a. Getting files ready for R
The goal is to make it into a phyloseq object that can be joined comptaible with our actual 18S rRNA dataset phyloseq object
```bash
conda activate phyloseq_tree
#this will be the taxonomy
cp ../reference_trees/annotations.koro.txt ~/legionella/R/phyloseq_tree/korotnevella
#metadata
sed -i '1 i\OTU apr13_E_73' ../reference_trees/annotations.koro.txt #this actually edits the file
#frequency data
sed 's/\t.*/\t 0/g' ../reference_trees/annotations.koro.txt  | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txt
mv ref_ASV.txt ~/legionella/R/phyloseq_tree/korotnevella
cd ~/legionella/R/phyloseq_tree/korotnevella
```
### 4b. Loading packages in R
```R
library(devtools)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ape)
library(microbiome)
library(ggpubr)
library(rstatix)
```
### 4c. Making phyloseq object for your ASVs
```R
#getting 18S ASVs
#reading sequence tables
seqtab_18S <- read.table("../../18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
#reading in taxonomy
tax_18S <- read.table("../../18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#reeading in metadata
map_18S <- read.table("../../metadata_18S.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
```
### 4d. Filtering your ASVs 
```R
physeq_18S_clr_filter = subset_taxa(physeq_18S, V7=="Korotnevella") #getting just specific genus of interest
ps1 <- merge_samples(physeq_18S, "month") #combine by month
physeq_18S_clr <- microbiome::transform(ps1, "clr") #converting abundance to CLR
```
### 4e. Getting your ASVs into a data frame
```R
data_18S <- psmelt(physeq_18S_clr_filter) #getting ASVs in datafrane
data_18S_2 <- data_18S[, c(1,2,3)] # selecting only columns of interest
```
### 4f. Reading in tree
```R
tre <- read.tree("korotnevella.tre")
```
### 4g. Get genome frequency table ready for phyloseq
```R
#reference otu table
seqtab_18S_1 <- read.table("../../18S-ASV-renamed.txt", header=T) # your ASV frequency table
seqtab_ref <- read.table("ref_ASV.txt", header=T) #reference genomes
#combing the otu tables
otu <- dplyr::full_join(seqtab_18S_1, seqtab_ref) #join reference and ASV frequency together (to get the same number of samples)
row.names(otu) <- otu$OTU
otu <- otu[,-1] # remove first row since it is a repeat of the row names now
otu <- data.matrix(otu) #convert to matrix
otu[is.na(otu)] <- as.numeric(0) #convert NAs to 0
otu_df <- as.data.frame(otu) #convert back into a dataframe
otu_df <- cbind(rownames(otu_df), data.frame(otu_df, row.names=NULL)) #make the rownames into the first column
colnames(otu_df)[1] <- "OTU" #rename the row names
filter_ref <- otu_df %>% 
  filter(!grepl('ASV', OTU)) #filter out all rows that start with ASV in first column
row.names(filter_ref) <- filter_ref$OTU #make first column into row names
filter_ref <- filter_ref[,-1] # remove first column
filter_ref <- data.matrix(filter_ref) #convert into matrix
otu_ref <-otu_table(filter_ref, taxa_are_rows=T) #make into phyloseq object
```
### 4h. Get other reference genome information into phyloseq object and merge
```R
#ref taxa
tax_other <- read.table("annotations.koro.txt", header=F, row.names=1, sep="\t")
tax_other_phylo <- tax_table(as.matrix(tax_other))
#metadata has already been read in
#combine ref data
physeq_ref <- merge_phyloseq(otu_ref,map_map_18S,tax_other_phylo,tre) #must have the tree added here or else funny NA thing happens in heatmap
physeq_ref_month <- merge_samples(physeq_ref, "month") #combine by month
```
### 4i. Get genome frequency into dataframe
```R
#getting ref genomes in df
data_ref <- psmelt(physeq_ref_month) #getting ref genomes into dataframe
data_ref_2 <- data_ref[, c(1,2,3)] #getting only important columns
data_ref_2[data_ref_2 == 0] <- NaN #replace 0 with NA
```
### 4j. Combine ref and ASV dataframes together
```R
all_genomes <- dplyr::full_join(data_18S_2, data_ref_2)
```
### 4k. Getting orders for plot
```R
month_order=c("march","april","may","june", "july","august") #month order
tip_order <- rev(tre$tip.label) #tip order
tip_order2 <- gsub("'","",tip_order) #cleaning up tip names
```
### 4l. Plot heatmap
```R
pdf("korotnevella_asv_heatmap_clr.pdf")
ggplot(all_genomes, aes(x=factor(Sample, levels=month_order),y=factor(OTU, levels=tip_order2))) +
  geom_tile((aes(fill = Abundance)))+
  scale_fill_gradient(low = "#FFC300", high = "#900C3F") +
  coord_fixed() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 2.5, vjust = 0.6, hjust = 0.0,), axis.text.y = element_text(size = 3),
        panel.background = element_rect(fill = 'grey67'), axis.title.x = element_text(size= 7), 
        axis.title.y = element_text(size= 6), plot.title = element_text(size = 10), 
        legend.title = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text = element_text(size = 4)) 
dev.off()
```
## 5. Month Boxplots
### 5a. Getting object ready for boxplot
```R
network_order=c("march_april","may_june","july_august") #order of bimonthly groups
month_colors=c("march"="darkorchid1","april"="#147BD1","may"="#2DC84D","june"="#F7EA48","july"="orange1","august"="#E03C31") #month colors
physeq_18S_bp_clr <- microbiome::transform(physeq_18S, "clr") #clr transformation
glom_bp <- tax_glom(physeq_18S_bp_clr, taxrank=rank_names(physeq_18S_bp_clr)[6]) #collapse by genus level
physeq_18S_clr_filter_bp = subset_taxa(glom_bp, V7=="Korotnevella") #get only genus of interest
```
### 5b. Converting phyloseq object into boxplot
```R
data_18S_bp <- psmelt(physeq_18S_clr_filter_bp) #getting object into a dataframe
data_pb_2 <- data_18S_bp[, c(2,29,3,6)] #getting only columns of interest (sample name, abundance, month, and bimonthly category)
```
### 5c. Plot boxplot
```R
pdf("koro_asv_boxplot_clr.pdf")
ggplot(data_pb_2, aes(x=factor(network, levels=network_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = month_colors)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
```
## 6. Combine tree and heatmap and boxplots
Use some photo editing software to align the tree and heatmap





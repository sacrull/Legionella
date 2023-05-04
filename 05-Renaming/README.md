# Renaming 18S ASVs
## 1. Load packages
```R
library(tidyverse)
library(reshape2)
library(phyloseq)
library(seqinr)
library(Biostrings)
```
## 2. Get into phyloseq object
Frequency table
```R
seqtab_18S <- read.table("feature-table-18S.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
```
Taxonomy
```R
tax_18S <- read.table("taxonomy_18S_clean.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
```
Representative Sequences
```R
refseq <- Biostrings::readDNAStringSet("18S_rep_seqs_silva.fasta")
```
Merge together all the phyloseq objects
```R
physeq_18S <- merge_phyloseq(otu_18S, tax_18S_phylo,refseq)
```
## 3. Rename
```R
taxa_names(physeq_18S) <- paste0("ASV", seq(ntaxa(physeq_18S)))
```
## 4. Get objects out of phyloseq and R
Get frequency table
```R
new_ASV <- t(as(otu_table(physeq_18S), "matrix"))
write.table(new_ASV, "18S-ASV.txt",col.names=NA,row.names=TRUE, sep= "\t") # says first column are rownames
```
Get taxonomy
```R
new_tax <- as(tax_table(physeq_18S), "matrix")
write.table(new_tax, "18S-tax.txt",col.names=FALSE,row.names=TRUE, sep= "\t") # says first column are rownames
```
Get representative sequences
```R
new_ref<-as.data.frame(refseq(physeq_18S))
refseq(physeq_18S) %>% Biostrings::writeXStringSet("plz.fna", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
```
## 5. Clean up new files in terminal
Fixing representative sequences
```bash
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' plz.fna > ref-seq-18S.fasta
```
Cleaning up taxa and frequency table
```bash
sed 's/"//g' 18S-tax.txt > 18S-tax-renamed.txt
sed 's/"//g' 18S-ASV.txt | sed -z 's/\t/OTU\t/'  > 18S-ASV-renamed.txt
```
## 6. Copy files to where you want them
```bash
cp 18S-tax-renamed.txt ~/legionella/R
cp 18S-ASV-renamed.txt ~/legionella/R
```

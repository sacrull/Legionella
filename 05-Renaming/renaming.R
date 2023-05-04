library(tidyverse)
library(reshape2)
library(phyloseq)
library(seqinr)
library(Biostrings)

seqtab_18S <- read.table("feature-table-18S.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
tax_18S <- read.table("taxonomy_18S_clean.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#rep sequences
refseq <- Biostrings::readDNAStringSet("18S_rep_seqs_silva.fasta")
#merge
physeq_18S <- merge_phyloseq(otu_18S, tax_18S_phylo,refseq)

#reanme
taxa_names(physeq_18S) <- paste0("ASV", seq(ntaxa(physeq_18S)))

#read out new objects
new_ASV <- t(as(otu_table(physeq_18S), "matrix"))
write.table(new_ASV, "18S-ASV.txt",col.names=NA,row.names=TRUE, sep= "\t") # says first column are rownames

new_tax <- as(tax_table(physeq_18S), "matrix")
write.table(new_tax, "18S-tax.txt",col.names=FALSE,row.names=TRUE, sep= "\t") # says first column are rownames

new_ref<-as.data.frame(refseq(physeq_18S))
refseq(physeq_18S) %>% Biostrings::writeXStringSet("plz.fna", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' plz.fna > ref-seq-18S.fasta

sed 's/"//g' 18S-tax.txt > 18S-tax-renamed.txt
sed 's/"//g' 18S-ASV.txt > 18S-ASV-renamed.txt
#manually add OTU to the first cell 18S-ASV-renamed.txt

cp 18S-tax-renamed.txt ~/legionella/R
cp 18S-tax-renamed.txt ~/legionella/eukref
cp 18S-tax-renamed.txt ~/legionella/eukref/vannella

cp 18S-ASV-renamed.txt ~/legionella/R
cp 18S-ASV-renamed.txt ~/legionella/eukref
cp 18S-ASV-renamed.txt ~/legionella/eukref/vannella
#getting databases
wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz
gzip -d *gz

#remove ward wrap for silva database
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' SILVA_138_SSURef_tax_silva.fasta > temp
mv temp SILVA_138_SSURef_tax_silva.fasta

#get vannella 
grep "Vannella" pr2_version_5.0.0_SSU_taxo_long.fasta -A 1 | sed 's/--//' > pr2_vannella.fa
grep "Vannella" SILVA_138_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_vannella.fa
ncbi query: Vannella 18S rRNA NOT whole genome shotgun NOT mrna 

sed '/^[^>]/s/U/T/g' silva_vannella.fa > temp
mv temp silva_vannella.fa

#remove word wrap for ncbi database and reformat header
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' ncbi_vannella.fa > temp
mv temp ncbi_vannella.fa

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' sequence_plat.fasta > temp
mv temp sequence_plat.fasta

#get all refrences together
cat ncbi_vannella.fa silva_vannella.fa pr2_vannella.fa sequence_plat.fasta > all_vannella.fa
#remove selectred genomes
mkdir remove
cp all_vannella.fa ./remove
cp remove.txt ./remove
sed 's/_.*//g' all_vannella.fa > temp
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' temp
for i in `cat remove.txt`; do rm ${i}.fasta; done
cat *.fasta > filt_vannella.fa
cp filt_vannella.fa ~/legionella/eukref/vannella/test1

#sort by and cluster by 99%
vsearch --sortbylength filt_vannella.fa --output all_vannella.sort
vsearch --cluster_fast all_vannella.sort --centroids all_vannella.clust --id 0.99

#clean up headers and use mafft
sed 's/|.*//' all_vannella.clust > temp
mv temp all_vannella.clust
mafft --auto all_vannella.clust > all_vannella.align.fa

#clean up for trimal
sed '/^[^>]/s/n/-/g' all_vannella.align.fa > temp
mv temp all_vannella.align.fa

#use trimal
~/bin/trimal -in all_vannella.align.fa -out all_vannella.trim.fa  -gt 0.3 -st 0.001

#make backbone tree
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n vannella_ref.tre -s all_vannella.trim.fa

#organize
mkdir reference_trees
grep ">" all_vannella.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.vannella.txt
sed 's/18S_rRNA_nucleus_strain_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' annotations.vannella.txt > temp1
sed 's/8S_rRNA_nucleus__Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp1 > temp2
sed 's/_18S_.*//g' temp2 > temp3
sed 's/18S_rRNA_nucleus_clone_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp3 > temp4
sed 's/_U//g' temp4 > annotations_vannella.txt
mv annotations_vannella.txt reference_trees/

mv RAxML_bipartitions.vannella_ref.tre reference_trees

awk '{print $1, $7}' ../18S-tax-renamed.txt | grep Vannella | awk '{print $1}' > vannella.ids

#pull their fasta file
cat vannella.ids | while read line; do grep -w $line ../ref-seq-18S.fasta -A 1; done > vannella.fa

#allign sequences
mkdir placement_tree
cd placement_tree

sina -i ../all_vannella.align.fa --prealigned -o all_vannella.arb
sina -i ../vannella.fa -r all_vannella.arb -o query.align.vannella.fa --fs-msc 0.01 --fs-full-len=100

#concat refrence sequences
cat query.align.vannella.fa ../all_vannella.align.fa > queryPlus.align.vannella.fa

#create tree
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n vannella.all.tre -s queryPlus.align.vannella.fa -t ../reference_trees/RAxML_bipartitions.vannella_ref.tre -T 7

#clean up tree
sed 's/QUERY___//g' RAxML_labelledTree.vannella.all.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.vannella.epa.tre

#get constraint tree
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n vannella.cons.tre -s queryPlus.align.vannella.fa -g ../reference_trees/RAxML_bipartitions.vannella_ref.tre -T 7 -x 25734 -p 25793

#export for midpoint rooting
RAxML_bipartitions.vannella.cons.tre

#make complete annotation file
awk '{print $1,$NF}' ../vannella.ids | sed 's/ /\t/g' | sed 's/ASV/ASV /2' > my_vannella_annots.txt
cat ../reference_trees/annotations_vannella.txt my_vannella_annots.txt > all_vannella_annotations.txt



#tricking phyloseq

cp ../reference_trees/annotations_vannella.txt ~/legionella/R/phyloseq_tree/vannella
sed -i '1 i\OTU apr13_E_73' ../reference_trees/annotations_vannella.txt #this actually edits the file
sed 's/\t.*/\t 0/g' ../reference_trees/annotations_vannella.txt  | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txt
mv ref_ASV.txt ~/legionella/R/phyloseq_tree/vannella

awk '{print $1, $7}' ../../18S-tax-renamed.txt > genus
cat genus ../reference_trees/annotations_vannella.txt  > taxa_all_vannella
sort taxa_all_vannella | uniq > taxa_vannella
sed 's/\t/\t /g' ../reference_trees/annotations_vannella.txt  > ref_taxa.txt

mv ref_ASV.txt ~/legionella/R/phyloseq_tree/vannella
#midpoint and order nodes RAxML_bipartitions.vannella.cons.tre then save and place in ~/legionella/R/phyloseq_tree/vannella/vannella.tre


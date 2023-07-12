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

#remove selectred genomes
mkdir remove
cd remove
sed 's/_.*//g' ../all_koro.fa > temp
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' temp
for i in `cat ../remove.txt`; do rm ${i}.fasta; done
cat *.fasta > filt_koro.fa
cd ..

#sort by and cluster by 99%
vsearch --sortbylength ./remove/filt_koro.fa --output all_koro.sort #when going back change to filt_acan.fa
vsearch --cluster_fast all_koro.sort --centroids all_koro.clust --id 0.99

#clean up headers and use mafft
sed 's/|.*//' all_koro.clust > temp
mv temp all_koro.clust
mafft --auto all_koro.clust > all_koro.align.fa

#clean up for trimal
sed '/^[^>]/s/n/-/g' all_koro.align.fa > temp
mv temp all_koro.align.fa

#use trimal
~/bin/trimal -in all_koro.align.fa -out all_koro.trim.fa  -gt 0.3 -st 0.001

#make backbone tree
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 7 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n koro_ref.tre -s all_koro.trim.fa

#organize
mkdir reference_trees
grep ">" all_koro.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.koro.txt
mv annotations.koro.txt reference_trees/
#sed 's/18S_rRNA_nucleus_strain_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' annotations.vannella.txt > temp1
#sed 's/8S_rRNA_nucleus__Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp1 > temp2
#sed 's/_18S_.*//g' temp2 > temp3
#sed 's/18S_rRNA_nucleus_clone_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp3 > temp4
#sed 's/_U//g' temp4 > annotations_vannella.txt

mv RAxML_bipartitions.koro_ref.tre reference_trees

awk '{print $1, $7}' ../18S-tax-renamed.txt | grep Korotnevella| awk '{print $1}' > koro.ids

#pull their fasta file
cat koro.ids | while read line; do grep -w $line ../ref-seq-18S.fasta -A 1; done > koro.fa

#allign sequences
mkdir placement_tree
cd placement_tree

sina -i ../all_koro.align.fa --prealigned -o all_koro.arb
sina -i ../koro.fa -r all_koro.arb -o query.align.koro.fa --fs-msc 0.01 --fs-full-len=100

#concat refrence sequences
cat query.align.koro.fa ../all_koro.align.fa > queryPlus.align.koro.fa

#create tree
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n koro.all.tre -s queryPlus.align.koro.fa -t ../reference_trees/RAxML_bipartitions.koro_ref.tre -T 7

#clean up tree
sed 's/QUERY___//g' RAxML_labelledTree.koro.all.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.koro.epa.tre

#get constraint tree
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n koro.cons.tre -s queryPlus.align.koro.fa -g ../reference_trees/RAxML_bipartitions.koro_ref.tre -T 7 -x 25734 -p 25793

#export for midpoint rooting
RAxML_bipartitions.koro.cons.tre

#make complete annotation file
awk '{print $1,$NF}' ../koro.ids | sed 's/ /\t/g' | sed 's/ASV/ASV /2' > my_koro_annots.txt
cat ../reference_trees/annotations.koro.txt my_koro_annots.txt > all_koro_annotations.txt



#tricking phyloseq

cp ../reference_trees/annotations.koro.txt ~/legionella/R/phyloseq_tree/korotnevella
sed -i '1 i\OTU apr13_E_73' ../reference_trees/annotations.koro.txt #this actually edits the file
sed 's/\t.*/\t 0/g' ../reference_trees/annotations.koro.txt  | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txt
mv ref_ASV.txt ~/legionella/R/phyloseq_tree/korotnevella

#midpoint and order nodes RAxML_bipartitions.vannella.cons.tre then save and place in ~/legionella/R/phyloseq_tree/vannella/vannella.tre

#awk '{print $1, $7}' ../../18S-tax-renamed.txt > genus
#cat genus ../reference_trees/annotations_vannella.txt  > taxa_all_vannella
#sort taxa_all_vannella | uniq > taxa_vannella
#sed 's/\t/\t /g' ../reference_trees/annotations_vannella.txt  > ref_taxa.txt










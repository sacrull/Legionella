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

#get acanthamoeba 
grep "Acanthamoeba" pr2_version_5.0.0_SSU_taxo_long.fasta -A 1 | sed 's/--//' > pr2_acan.fa
grep "Acanthamoeba" SILVA_138_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_acan.fa
ncbi query: Acanthamoeba 18S rRNA NOT whole genome shotgun NOT mrna 

sed '/^[^>]/s/U/T/g' silva_acan.fa > temp
mv temp silva_acan.fa

#remove word wrap for ncbi database and reformat header
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' sequence_acan.fasta > temp
mv temp ncbi_acan.fa

#get all refrences together
cat ncbi_acan.fa silva_acan.fa pr2_acan.fa > all_acan.fa



#remove selectred genomes
mkdir remove
cp all_acan.fa ./remove
cp remove.txt ./remove
cd remove
sed 's/_.*//g' all_acan.fa > temp
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' temp
for i in `cat remove.txt`; do rm ${i}.fasta; done
cat *.fasta > filt_acan.fa
cp filt_acan.fa ~/legionella/eukref/acanthamoeba/
cd ..

#sort by and cluster by 99%
vsearch --sortbylength ./remove/filt_acan.fa --output all_acan.sort #when going back change to filt_acan.fa
vsearch --cluster_fast all_acan.sort --centroids all_acan.clust --id 0.99

#clean up headers and use mafft
sed 's/|.*//' all_acan.clust > temp
mv temp all_acan.clust
mafft --auto all_acan.clust > all_acan.align.fa

#clean up for trimal
sed '/^[^>]/s/n/-/g' all_acan.align.fa > temp
mv temp all_acan.align.fa

#use trimal
~/bin/trimal -in all_acan.align.fa -out all_acan.trim.fa  -gt 0.3 -st 0.001

#make backbone tree
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n acan_ref.tre -s all_acan.trim.fa

#organize
mkdir reference_trees
grep ">" all_acan.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.acan.txt
mv annotations.acan.txt reference_trees/
#sed 's/18S_rRNA_nucleus_strain_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' annotations.vannella.txt > temp1
#sed 's/8S_rRNA_nucleus__Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp1 > temp2
#sed 's/_18S_.*//g' temp2 > temp3
#sed 's/18S_rRNA_nucleus_clone_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp3 > temp4
#sed 's/_U//g' temp4 > annotations_vannella.txt

mv RAxML_bipartitions.acan_ref.tre ./reference_trees

#skipping unassigned ASVs since only interested in the assigned ones
#create a list of all ASV ids vannella_ASV_names
awk '{print $1, $7}' ../18S-tax-renamed.txt | grep Acanthamoeba| awk '{print $1}' > acan.ids

#pull their fasta file
cat acan.ids | while read line; do grep -w $line ../ref-seq-18S.fasta -A 1; done > acan.fa

#allign sequences
mkdir placement_tree
cd placement_tree

sina -i ../all_acan.align.fa --prealigned -o all_acan.arb
sina -i ../acan.fa -r all_acan.arb -o query.align.acan.fa --fs-msc 0.01 --fs-full-len=100

#concat refrence sequences
cat query.align.acan.fa ../all_acan.align.fa > queryPlus.align.acan.fa

#create tree
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n acan.all.tre -s queryPlus.align.acan.fa -t ../reference_trees/RAxML_bipartitions.acan_ref.tre -T 7

#clean up tree
sed 's/QUERY___//g' RAxML_labelledTree.acan.all.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.acan.epa.tre

#get constraint tree
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n acan.cons.tre -s queryPlus.align.acan.fa -g ../reference_trees/RAxML_bipartitions.acan_ref.tre -T 7 -x 25734 -p 25793

#make complete annotation file
awk '{print $1,$NF}' ../acan.ids | sed 's/ /\t/g' | sed 's/ASV/ASV /2' > my_acan_annots.txt
cat ../reference_trees/annotations.acan.txt my_acan_annots.txt > all_acan_annotations.txt

#export for midpoint rooting
RAxML_bipartitions.acan.cons.tre

cp ../reference_trees/annotations.acan.txt ~/legionella/R/phyloseq_tree/acanthamoeba
sed -i '1 i\OTU apr13_E_73' ../reference_trees/annotations.acan.txt #this actually edits the file
sed 's/\t.*/\t 0/g' ../reference_trees/annotations.acan.txt  | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txtf
mv ref_ASV.txt ~/legionella/R/phyloseq_tree/acanthamoeba
#midpoint and order nodes RAxML_bipartitions.acan.cons.tre then save and place in ~/legionella/R/phyloseq_tree/acanthamoeba/acan.tre










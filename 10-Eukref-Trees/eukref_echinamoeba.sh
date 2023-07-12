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

#get Echinamoebida 
grep "Echinamoebida" pr2_version_5.0.0_SSU_taxo_long.fasta -A 1 | sed 's/--//' > pr2_echi.fa
grep "Echinamoebida" SILVA_138_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_echi.fa
ncbi query: Echinamoebida 18S rRNA NOT whole genome shotgun NOT mrna 

sed '/^[^>]/s/U/T/g' silva_echi.fa > temp
mv temp silva_echi.fa

#remove word wrap for ncbi database and reformat header
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' sequenc_echi.fasta > temp
mv temp ncbi_echi.fa

#get all refrences together
cat ncbi_echi.fa silva_echi.fa pr2_echi.fa > all_echi.fa

#remove selectred genomes
mkdir remove
grep ">" all_echi.fa | grep -v Echinamoeba | sed 's/ .*//g' | sed 's/>//' > remove.txt
cd remove
sed 's/_.*//g' ../all_echi.fa > temp
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' temp
for i in `cat ../remove.txt`; do rm ${i}.fasta; done
cat *.fasta > filt_echi.fa
cd ..

#sort by and cluster by 99%
vsearch --sortbylength ./remove/filt_echi.fa --output all_echi.sort #when going back change to filt_acan.fa
vsearch --cluster_fast all_echi.sort --centroids all_echi.clust --id 0.99

#clean up headers and use mafft
sed 's/|.*//' all_echi.clust > temp
mv temp all_echi.clust
mafft --auto all_echi.clust > all_echi.align.fa

#clean up for trimal
sed '/^[^>]/s/n/-/g' all_echi.align.fa > temp
mv temp all_echi.align.fa

#use trimal
~/bin/trimal -in all_echi.align.fa -out all_echi.trim.fa  -gt 0.3 -st 0.001

#make backbone tree
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 7 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n echi_ref.tre -s all_echi.trim.fa

#organize
mkdir reference_trees
grep ">" all_echi.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' | sed 's/_U//' > annotations.echi.txt
mv annotations.echi.txt reference_trees/
#grep ">" all_echi.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' | grep Echinamoebida | grep -v vermiformis > annotations.echi.txt
#sed 's/18S_rRNA_nucleus_strain_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' annotations.vannella.txt > temp1
#sed 's/8S_rRNA_nucleus__Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp1 > temp2
#sed 's/_18S_.*//g' temp2 > temp3
#sed 's/18S_rRNA_nucleus_clone_.*_Eukaryota_Amoebozoa_Discosea_Discosea_X_Flabellinia_Vannellida_Vannellidae_//g' temp3 > temp4
#sed 's/_U//g' temp4 > annotations_vannella.txt

mv RAxML_bipartitions.echi_ref.tre reference_trees

awk '{print $1, $7}' ../18S-tax-renamed.txt | grep Echinamoebida | awk '{print $1}' > echi.ids

#pull their fasta file
cat echi.ids | while read line; do grep -w $line ../ref-seq-18S.fasta -A 1; done > echi.fa


#allign sequences
mkdir placement_tree
cd placement_tree

sina -i ../all_echi.align.fa --prealigned -o all_echi.arb
sina -i ../echi.fa -r all_echi.arb -o query.align.echi.fa --fs-msc 0.01 --fs-full-len=100

#concat refrence sequences
cat query.align.echi.fa ../all_echi.align.fa > queryPlus.align.echi.fa

#create tree
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n echi.all.tre -s queryPlus.align.echi.fa -t ../reference_trees/RAxML_bipartitions.echi_ref.tre -T 7

#clean up tree
sed 's/QUERY___//g' RAxML_labelledTree.echi.all.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.echi.epa.tre

#get constraint tree
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n echi.cons.tre -s queryPlus.align.echi.fa -g ../reference_trees/RAxML_bipartitions.echi_ref.tre -T 7 -x 25734 -p 25793

#export for midpoint rooting
RAxML_bipartitions.echi.cons.tre

#make complete annotation file
awk '{print $1,$NF}' ../echi.ids | sed 's/ /\t/g' | sed 's/ASV/ASV /2' > my_echi_annots.txt
cat ../annotations.echi.txt my_echi_annots.txt > all_echi_annotations.txt

#tricking phyloseq

cp ../reference_trees/annotations.echi.txt ~/legionella/R/phyloseq_tree/echinamoebida
sed -i '1 i\OTU apr13_E_73' ../reference_trees/annotations.echi.txt #this actually edits the file
sed 's/\t.*/\t 0/g' ../reference_trees/annotations.echi.txt  | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txt
mv ref_ASV.txt ~/legionella/R/phyloseq_tree/echinamoebida

#midpoint and order nodes RAxML_bipartitions.echi.cons.tre then save and place in ~/legionella/R/phyloseq_tree/echinamoebida/echi.tre

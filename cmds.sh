# make a directory for hw5

mkdir hw5 

cd hw5 

#use conda environment 

conda activate hw5

#download the files for the hw

fasterq-dump \
 SRR3214715 SRR3215024 SRR3215107 

# install fastp if not there, I have already installed fastp and other packages for the work that needs to be done in this 


fastp \
 -i SRR3214715_1.fastq \
 -I SRR3214715_1_trimmed.fastq
 -O SRR3214715_2_trimmed.fastq \
 --html SRR3214715_fastp_report.html

# same can steps needs to be doen for all the files, e.g SRR3215024 SRR3215107



skesa \
 --reads SRR3214715_1_trimmed.fastq SRR3214715_2_trimmed.fastq \
 --contigs_out SRR3214715_assembly.fna

# this step needs to be done again for all the trimmed files

#filter contigs using the .py script provided

python ./filter.contigs.py \
 -i SRR3214715_assembly.fna \
 -l 1000 \
 -c 5 \
 -m \
 -o SRR3214715_filtered_assembly.fna 

# if the script is not working, do chmod+x ./filter_contigs.py

# this step also needs to be done for all the files


#check all assemblies are approximately the same size
ls -alh *.fna

#make a directory mlst

mkdir mlst


# run MLST with docker
docker pull staphb/mlst:latest
docker run -it --mount type=bind,src=$HOME/hw5,target=/local staphb/mlst bash
cd /local
mlst SRR3214715_filtered_assembly.fna > ./mlst/SRR3214715_MLST_Summary.tsv
mlst SRR3215024_filtered_assembly.fna > ./mlst/SRR3215024_MLST_Summary.tsv
mlst SRR3215107_filtered_assembly.fna > ./mlst/SRR3215107_MLST_Summary.tsv
exit

#concatenate all the summary files

cat /hw5/mlst/*_MLST_Summary.tsv > mlst_noheaders.tsv


# just get the first column with accession id: 

awk -F'[ /]' '{print $2, $4, $5, $6, $7, $8, $9, $10}' mlst_noheaders.tsv | sed 's/_filtered_assembly.fna//g' | sed 's/campylobacter//g' > mlst.tsv


#make a directory and move to fastani
mkdir fastani 
cd fastani

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/821/535/GCF_017821535.1_ASM1782153v1/GCF_017821535.1_ASM1782153v1_genomic.fna.gz
gunzip *.fna
mv GCF_017821535.1_ASM1782153v1_genomic.fna reference.fna

# run the command for fastani

fastANI \
  --query /hw5/SRR3214715_filtered_assembly.fna \
  --ref reference.fna \
  --fragLen 1000 \
  --output SRR3214715_fastANI_Output.tsv

# do this step for all the filtered files 


# combine tsv files
cat *_fastANI_Output.tsv > FastANI_Output.tsv


#calculate alignment % and len

awk \
  '{alignment_percent = $4/$5*100} \
   {alignment_length = $4*3000} \
   {print $0 "\t" alignment_percent "\t" alignment_length}' \
  FastANI_Output.tsv \
  > FastANI_Output_With_Alignment.tsv
sed \
  "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" \
  FastANI_Output_With_Alignment.tsv \
  > FastANI_Output_With_Alignment_With_Header.tsv
column -ts $'\t' FastANI_Output_With_Alignment_With_Header.tsv | less -S


# replace query names with accession IDs
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0} NR>1 {match($1, /SRR[0-9]+/); print substr($1, RSTART, RLENGTH), $2, $3, $4, $5, $6, $7}' FastANI_Output_With_Alignment_With_Header.tsv > FastANI_Output_With_Alignment_With_Header_modified.tsv


# make a directory and change it

mkdir checkm

cd checkm

#make a directory asm and put the any one of the assembled files in this
mkdir asm
mv /hw5/SRR3214715_assembly.fna /hw5/cheeckm


# get the database: 

# get checkm database
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz
tar zxvf checkm_data_2015_01_16.tar.gz
echo 'export CHECKM_DATA_PATH=/home/amckee/biol7210/ex5/checkm/db' >> ~/.bashrc
source ~/.bashrc
echo "${CHECKM_DATA_PATH}"

checkm taxon_list | grep Campylobacter
checkm taxon_set species "Campylobacter jejuni" Cj.markers
checkm \
  analyze \
  Cj.markers \
  ~/hw5/checkm/asm \
  analyze_output
checkm \
  qa \
  -f checkm.tax.qa.out \
  -o 1 \
  Cj.markers \
  analyze_output
sed 's/ \+ /\t/g' checkm.tax.qa.out > checkm.tax.qa.out.tsv
cut -f 2- checkm.tax.qa.out.tsv > tmp.tab && mv tmp.tab checkm.tax.qa.out.tsv
sed -i '1d; 3d; $d' checkm.tax.qa.out.tsv

# extract the first 
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0} NR>1 {sub(/_filtered_assembly$/, "", $1); print}' checkm.tax.qa.out.tsv > quality.tsv




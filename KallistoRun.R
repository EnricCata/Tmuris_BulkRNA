# -*- coding: utf-8 -*-
# trichuris_muris_genomics_Life_Cycle.Bulkrnaseq

# ## Wenxin Jiang & Enric Cata

# - code used to explore bulk RNAseq data from Maria Duque's lab investigating gene expression/regu;ation key in T. muris development(molting) in mice and in organoids



## Download reference from Sanger
```bash
# working directory
cd /home/wenxinjiang/Lab_data/Bulk_RNAseq_Life_Cycle/T.muris_reference

# genome
wget https://ngs.sanger.ac.uk/production/pathogens/sd21/reference_genomes/trichuris_muris/trichuris_muris_wsi-v7.0.fa

# annotation
wget https://ngs.sanger.ac.uk/production/pathogens/sd21/reference_genomes/trichuris_muris/trichuris_muris.annotation.v7.1.gff3.gz  
```



## Get the raw RNAseq data
- Stephen Doyle (Steve) share me with metadata of samples all from the lane 46825_1
- shared under link: ftp://ngs.sanger.ac.uk/production/pathogens/sd21/tmuris_rnaseq
- need to download these into my server, they are fastq format
- A lanes_samples.list was also shared to me which connects the file names to the sample names.

```bash

# download data in my server from shared directory (%23 = # as system does not recognize #; following are the example commend for #36 file)
wget ftp://ngs.sanger.ac.uk/production/pathogens/sd21/tmuris_rnaseq/46825_1%2336.fastq.gz
wget ftp://ngs.sanger.ac.uk/production/pathogens/sd21/tmuris_rnaseq/46825_1%2336_1.fastq.gz
wget ftp://ngs.sanger.ac.uk/production/pathogens/sd21/tmuris_rnaseq/46825_1%2336_2.fastq.gz

# do the similar commend for each interested samples
```



# # prior to running callisto, TRIM!

# # Run Kallisto

```bash 
# make a transcripts fasta
# download gffread-compare module
gffread -x TRANSCRIPTS.fa -g REF.fa ANNOTATION.gff3
#结果：FASTA index file REF.fa.fai created.

# index the transcripts
./kallisto index --index TRANSCRIPTS.ixd TRANSCRIPTS.fa



# run kallisto for trimmed
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
kallisto quant \
  --bias \
  --index /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/06_scripts/kallisto/TMURIS-INDEX.idx \
  --output-dir /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/kallisto_mapped_samples/${SAMPLE}_kallisto_out \
  --bootstrap-samples 100 \
  --threads 7 \
  --fusion \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-trimmed/${SAMPLE}_1_val_1.fq.gz \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-trimmed/${SAMPLE}_2_val_2.fq.gz;
done < /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/SAMPLE_ID_NONADULT.txt


# read LANE SAMPLE_ID (this is a list or txt that contains the name, which gets invoked in the {LANE})
# test for one file
# check what each part does IN THE MANUAL!


# kallisto for non trimmed   
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
kallisto quant \
  --bias \
  --index /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/06_scripts/kallisto/TMURIS-INDEX.idx \
  --output-dir /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/kallisto_mapped_samples/${SAMPLE}_kallisto_out \
  --bootstrap-samples 100 \
  --threads 7 \
  --fusion \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/${SAMPLE}_1.fastq.gz  \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/${SAMPLE}_2.fastq.gz;
done < /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/SAMPLE_ID_ADULT.txt


kallisto quant \
  --bias \
  --index /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/06_scripts/kallisto/TMURIS-INDEX.idx \
  --output-dir /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/kallisto_mapped_samples/49707_1#23_kallisto_out \
  --bootstrap-samples 100 \
  --threads 7 \
  --fusion \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/49707_1#23_1.fastq.gz  \
  /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/02_raw-data/fastq-raws/49707_1#23_2.fastq.gz;



mkdir KALLISTO_MAPPED_SAMPLES
mv kallisto_* KALLISTO_MAPPED_SAMPLES/

#结果：output files: abundance.h5  abundance.tsv  fusion.txt  run_info.json
```


```bash

# extract TPMs per sample
# raw counts: better if you compare results with single RNA data
# TPM (transcript per million, normalizing transcript) 
for i in ` ls -1d *out `; do 
  echo $i > ${i}.tpm ; cat ${i}/abundance.tsv | cut -f5 | sed '1d' >> ${i}.tpm; 
  done
# cut and sed are removing junk reads from the file, additional info not relevant to the data (remove the 'labels')



# generate a "transcripts" list, taken from the TRANSCRIPTS.fa file
# echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d  " " | sed 's/>//g' >> transcripts.list
# due to Apollo giving long unique codes, the transcript IDs are obscure. Here is the fix
# awk '$3=="mRNA" {print $9}' ../ANNOTATION.gff3 | cut -f3,5 -d";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDtoNAME_conversion.txt

# while read ID NAME; do sed -i "s/${ID}/${NAME}/g" transcripts.list; done < mRNA_IDtoNAME_conversion.txt &

# ALTERNATE WAY, direct from the annotaiton

echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d" " | sed -e 's/>//g' >> transcripts.list

echo "ID" > transcripts.list; grep ">" /data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/06_scripts/kallisto/TMURIS-TRANSCRIPTS.fa | cut -f1 -d" " | sed -e 's/>//g' >> transcripts.list

# cuts and pastes everyhtng BUT the first line, the sed part retains punctuation (like dots, slashes, and question marks)

# make a data frame containing all TMP values from all samples
cp transcripts.list tmp; ls -1v *_kallisto_out.tpm | while read SAMPLE; do paste tmp $SAMPLE > tmp2; mv tmp2 tmp; done
# above is merging all TPM values extracted from each sample (*_out.tpm). Since its stored in separate files for each sample, merge them all together with one column per sample and ID for rows for everything!

mv tmp kallisto_allsamples.tpm.table

sed -i -e 's/Transcript://g' -e 's/kallisto_//g' -e 's/_out//g' kallisto_allsamples.tpm.table

# Wenxin used the second one since the first way doesn't quite work

```


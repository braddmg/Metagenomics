# Metagenomics
Pipeline for metagnomic analyses of water samples

## Quality control
First step is to remove wuality control reads. For this we used fastp.
Asumming all your samples (in fastq.gz format) are in the same folder.
```bash
for i in `ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//'`; do  fastp -i $i\_1.fastq.gz -I $i\_2.fastq.gz --detect_adapter_for_pe -o trimmed/$i\_1.fq.gz -O trimmed/$i\_2.fq.gz -h trimmed/$i\_fastq.html -e 25
```
New data will be saved in a new folder named trimmed.

Then use bowtie2 to remove reads associated with the human genome. Please review this link to download the indexed human genome. https://benlangmead.github.io/aws-indexes/bowtie
```bash
for i in `ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//'`; do
bowtie2 -p 8 -x databases/references/H.sapiens_hg19/bowtie2/hg19 \
  -1 trimmed/$i\_1.fq.gz \
  -2 trimmed/$i\_2.fq.gz \
--very-sensitive-local \
--quiet \
  --un-conc-gz filtered/$i\_filtered
done
```
Samples will be saved in a new folder named filtered.

## Taxonomic assignment and count of reads
For this purpose we will use kraken2 to asign taxonomy to the reads of each sample. 
Move into the filtered folder. 
Please see this link to select and download an specific database to assign taxonomy. In our case we used a custom database. https://benlangmead.github.io/aws-indexes/k2
```bash
for i in `ls -1 *_filtered.1.fq.gz | sed 's/_filtered.1.fq.gz//' `; do
kraken2 --paired $i\_filtered.1.fq.gz $i\_filtered.2.fq.gz --classified-out $i\_gen#.fq --report $i\_gen.kreport --db databases/gtdb/2023-04-25/databases/gtdb_r207_v2_genomes/kraken2/gtdb_r207_v2_genomes/ --threads 36 --gzip-compressed --confidence 0.5; done
```
Kraken2 generates reports in kreport extensio. We will use this files to create a matrix of read counts for each taxon across the samples wit bracken. 
```bash
for i in `ls -1 *.kreport | sed 's/.kreport//' `; do bracken -r 100 -i $i\.kreport -o bracken/$i\.bracken -d /data/databases/gtdb/2023-04-25/databases/gtdb_r207_v2_genomes/bracken; done
```
In the folder named bracken, we will have a .bracken file for ech sample. Those reports can be merged into a single file with the script combine_bracken_outputs: https://github.com/jenniferlu717/Bracken/blob/master/analysis_scripts/combine_bracken_outputs.pyhttps://github.com/jenniferlu717/Bracken/blob/master/analysis_scripts/combine_bracken_outputs.py

```bash
python combine_bracken_outputs.py --files *.bracken -o bracken_results.tsv
```
The result of this file will be a matrix where rowas are taxa associated to each read, the columns are the samples and the values is the total and fraction of reads assigned. 
We can convert this file into a phyloseq object in R to performe diversity analyses. 
Please see https://github.com/braddmg/b2p for more information. And additional microbial diversity analyses can be found here: 

## Metagenomic Coassembly

filtered samples will be used to generate a coassembly/assembly of metagenomes. 
We can use Megahit for this purpose and here is an example of a coassembly with three different samples, using both forward and reverse files. 
```bash
megahit -1 S035_filtered.1.fq.gz,S036_filtered.1.fq.gz,S037_filtered.1.fq.gz -2  S035_filtered.2.fq.gz,S036_filtered.2.fq.gz,S037_filtered.2.fq.gz --k-list 33,55,77,99,127 -m 512000000000 -t 32 -o Megahit/S1
```
The result (contigs) will be saved in the folder Megahit/S1

## Evaluating coassembly

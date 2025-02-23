# Metagenomics
Pipeline for metagenomic analyses of water samples

## Quality control
The first step in metagenomic data analysis is to perform quality control on the raw reads. For this, we use [fastp](https://github.com/OpenGene/fastp), an efficient tool for quality filtering, trimming adapters, and generating quality reports.
Asumming all your samples (in fastq.gz format) are in the same folder.
```bash
for i in `ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//'`; do  fastp -i $i\_1.fastq.gz -I $i\_2.fastq.gz --detect_adapter_for_pe -o trimmed/$i\_1.fq.gz -O trimmed/$i\_2.fq.gz -h trimmed/$i\_fastq.html -e 25
```
New data will be saved in a new folder named trimmed.

Then use bowtie2 to remove reads associated with the human genome. Please review this [link](https://benlangmead.github.io/aws-indexes/bowtie)to download the indexed human genome. 
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
Please see this [link](https://benlangmead.github.io/aws-indexes/k2) to select and download an specific database to assign taxonomy. In our case we used a custom database previously indexed with genomes from the [GTDB](https://gtdb.ecogenomic.org). 
```bash
for i in `ls -1 *_filtered.1.fq.gz | sed 's/_filtered.1.fq.gz//' `; do
kraken2 --paired $i\_filtered.1.fq.gz $i\_filtered.2.fq.gz --classified-out $i\_gen#.fq --report $i\_gen.kreport --db databases/gtdb/2023-04-25/databases/gtdb_r207_v2_genomes/kraken2/gtdb_r207_v2_genomes/ --threads 36 --gzip-compressed --confidence 0.5; done
```
Kraken2 generates reports in kreport extension. We will use this files to create a matrix of read counts for each taxon across the samples with bracken. 
```bash
for i in `ls -1 *.kreport | sed 's/.kreport//' `; do bracken -r 100 -i $i\.kreport -o bracken/$i\.bracken -d /data/databases/gtdb/2023-04-25/databases/gtdb_r207_v2_genomes/bracken; done
```
In the folder named bracken, we will have a .bracken file for ech sample. Those reports can be merged into a single file with the script [combine_bracken_outputs.py](https://github.com/jenniferlu717/Bracken/tree/master/analysis_scripts) 

```bash
python combine_bracken_outputs.py --files *.bracken -o bracken_results.tsv
```
The result of this file will be a matrix where rowas are taxa associated to each read, the columns are the samples and the values is the total and fraction of reads assigned. 
We can convert this file into a phyloseq object in R to performe diversity analyses. 
Please see [b2p github repository](https://github.com/braddmg/b2p) for more information. 

## Metagenomic Coassembly

Filtered samples will be used to generate a coassembly/assembly of metagenomes. 
We can use Megahit for this purpose and here is an example of a coassembly with three different samples, using both forward and reverse files. 
```bash
megahit -1 S035_filtered.1.fq.gz,S036_filtered.1.fq.gz,S037_filtered.1.fq.gz -2  S035_filtered.2.fq.gz,S036_filtered.2.fq.gz,S037_filtered.2.fq.gz --k-list 33,55,77,99,127 -m 512000000000 -t 32 -o Megahit/S1
```
The result (contigs) will be saved in the folder Megahit/S1

## Evaluating coassembly

Assembled contigs can be evaluated using MetaQuast, which provides quality assessment based on alignment and assembly metrics. 
To process all samples efficiently, locate the FASTA files (.fasta) and use the following command within a loop:
```bash
for i in `ls -1 *.fasta | sed 's/.fasta//'`
do
metaquast.py -L -s $i\.fasta -o QUAST/$i/ --min-contig 500
done
```
## Identifiying plasmidic contigs

There are multiple tools available for identifying plasmidic contigs in metagenomic data. In this workflow, we will use PlasX, a machine learning-based software designed for plasmid detection:
[PlasX GitHub Repository](https://github.com/michaelkyu/PlasX)
The following command processes all samples in FASTA format, filters out contigs shorter than 500 bp, annotates them using Anvi’o with the COG14 and Pfam databases, and assigns a PlasX score to each contig: 
```bash
for i in `ls *fasta | awk 'BEGIN{FS=".fasta"}{print $1}'`
do
anvi-script-reformat-fasta $i.fasta \
                           -o $i.fa \
                           -l 500 --seq-type NT --simplify-names --prefix $i
anvi-gen-contigs-database -f $i.fa -o $i.db
anvi-run-hmms -c $i.db
anvi-export-gene-calls --gene-caller prodigal -c $i.db -o $i-gene-calls.txt
done

for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
anvi-run-ncbi-cogs -T 32 --cog-version COG14 --cog-data-dir /work/databases/anvio/COG_2014 -c $i.db
anvi-run-pfams -T 32 --pfam-data-dir /work/bmendoza/Tesis/Data/plasmids/anvio/Pfam_v32 -c $i.db
anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c $i.db -o $i-cogs-and-pfams.txt
done

for i in `ls *fasta | awk 'BEGIN{FS=".fasta"}{print $1}'`
do
plasx search_de_novo_families \
    -g $i-gene-calls.txt \
    -o $i-de-novo-families.txt \
    --threads 32 \
    --splits 32 \
    --overwrite

plasx predict \
    -a $i-cogs-and-pfams.txt $i-de-novo-families.txt \
    -g $i-gene-calls.txt \
    -o $i-scores.txt \
    --overwrite
done
```
## Detecting antibiotic resistance genes (ARGs)
We can detect ARGs from fasta files with [ABRicate](https://github.com/tseemann/abricate)
In this example, we employ the [CARD](https://card.mcmaster.ca) database, ABRicate supports several other databases. 

```bash
for i in `ls -1 *.fa | sed 's/.fa//'`
do
abricate $i\.fa --db card > $i\.card
```
## Quatifying normalized abundance of ARGs in co-assemblies
First, index the metagenomic contigs and map the sample reads against these indexed contigs.

```bash
#Create a mapping file of each metagnenome.
bowtie2-build S1_M1.fa mapping/S1_M1
#Mapping each sample reads agaisnt the metagenome.
bowtie2 --threads 32 -x mapping/S1_M1 -1 ../QC_samples/MS22_filtered.1.fq.gz -2 ../QC_samples/MS22_filtered.2.fq.gz -S mapping/S1_M1_1.sam
bowtie2 --threads 32 -x mapping/S1_M1 -1 ../QC_samples/MS29_filtered.1.fq.gz -2 ../QC_samples/MS29_filtered.2.fq.gz -S mapping/S1_M1_2.sam
bowtie2 --threads 32 -x mapping/S1_M1 -1 ../QC_samples/MS36_filtered.1.fq.gz -2 ../QC_samples/MS36_filtered.2.fq.gz -S mapping/S1_M1_3.sam
```
We can now convert sam to bam files.
```bash
for i in `ls -1 ../QC_samples/*.sam | sed 's/.sam//'`
do
samtools view -F 4 -bS mapping/$i\.sam > mapping/$i\-RAW.bam
anvi-init-bam mapping/$i\-RAW.bam -o mapping/$i\.bam
rm mapping/$i\.sam mapping/$i\-RAW.bam
done
```
Use the R script below to extract the coordinates of ARG-containing regions from a CSV file (resistance.csv) annotated by ABRicate. For each sample, a region file (*_regions.txt) is created.
```R
#Load libraries
library(dplyr)
library(readr)
# Input file
csv_file <- "resistencia.csv"
# Read the CSV file
data <- read_csv(csv_file, col_types = cols(
  ID = col_character(),
  SEQUENCE = col_character(),
  START = col_integer(),
  END = col_integer()
))

# Loop through unique IDs and create separate regions.txt files
unique_ids <- unique(data$ID)
for (id in unique_ids) {
  # Filter rows for the current ID
  subset_data <- data %>% filter(ID == id)
  # Select and rename columns to match regions.txt format
  regions <- subset_data %>% 
    select(SEQUENCE, START, END) %>%
    rename(Contig = SEQUENCE, Start = START, End = END)
  # Write to regions.txt file
  output_file <- paste0(id, "_regions.txt")
  write_delim(regions, output_file, delim = "\t", col_names = FALSE)
  cat("Created:", output_file, "\n")
}
```
Use the following Python script to count reads aligned to ARG regions and compute RPKM values. You will need [pysam](https://github.com/pysam-developers/pysam) and [samtools](https://github.com/samtools/samtools). 
```Python
for i in `ls -1 *.bam | sed 's/.bam//'`
do
  python3 - <<EOF
import pysam
i = "$i"
regions_file = "${i}_regions.txt"
bam_file = "${i}.bam"
bam = pysam.AlignmentFile(bam_file, "rb")

# Function to count reads in a specific region
def count_reads(bam, contig, start, end):
    count = 0
    for read in bam.fetch(contig, start, end):
        if not read.is_unmapped:  # Exclude unmapped reads
            count += 1
    return count
# Calculate total mapped reads in BAM file
def total_mapped_reads(bam):
    return sum(1 for read in bam if not read.is_unmapped)
# Get total mapped reads
total_reads = total_mapped_reads(bam)

# Read regions from file and save results
output_file = f"{i}_read_counts_with_rpkm.txt"
with open(regions_file, "r") as regions, open(output_file, "w") as output:
    output.write("Contig\tStart\tEnd\tRead_Count\tRPKM\tTotal_Mapped_Reads\n")
    for line in regions:
        contig, start, end = line.strip().split("\t")
        start, end = int(start), int(end)
        region_length = end - start + 1
        # Count reads for the region
        read_count = count_reads(bam, contig, start, end)
        # Calculate RPKM
        if total_reads > 0:
            rpkm = (read_count * 1e9) / (region_length * total_reads)
        else:
            rpkm = 0
        # Save results to file
        output.write(f"{contig}\t{start}\t{end}\t{read_count}\t{rpkm:.6f}\t{total_reads}\n")

# Close BAM file
bam.close()
print(f"Processed {bam_file}, results saved to {output_file}")
EOF
done
```
## Assigning taxonomy to metagenomic contigs
For assigning taxonomy to metagenomic contigs we will use [Kaiju](https://github.com/bioinformatics-centre/kaiju) tool. This tool assigns the closest taxonomy to each contig in a metagenomic assembly using protein-level classification and various databases. 
In this example we will employ the SwissProt database. For instructions on indexing or downloading databases, refer to the [Kaiju github repository](https://github.com/bioinformatics-centre/kaiju).
```bash
kaiju -z 32 -t database/swissprot_nodes.dmp -f database/kaiju_db_refseq.fmi -i metagenomes.fa -o kaiju.out -v
```




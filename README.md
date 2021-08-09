# Nextflow Next-generation Sequence-analysis Toolkit (Nf-NeST) : A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data

This Nf-NeST is a nextflow-Docker version of Next-generation Sequence-analysis Toolkit(NeST) with improved snpfilter function
: A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data

1. [Overview of NeST framework](#Overview)
2. [Prerequisites](#Prerequisites)
3. [Availability of code and installation](#Installation)
4. [Your first analysis](#First)
5. [Input standardization](#inputs)
6. [Output Description](#outputs)


<a id="Overview"></a>
## Overview of the NeST framework:

NeST is a python based modular framework for consensus based variant calling. The overall analysis framework is broken down into four major blocks.
1. PrepInputs
2. Align reads and variant calling
3. SNPfiltering
4. Summarize report with json/PowerBI input

<a id="Prerequisites"></a>

## Prerequisites
- MAC
    - Download Docker Desktop here https://www.docker.com/get-started
- Ubuntu
    - https://docs.docker.com/engine/install/ubuntu/
    
<a id="Installation"></a>
## Availability of code and installation:

1. Download git repository:

   Clone the master branch of this repository.
   ```
   git clone https://github.com/CDCgov/Nf-NeST.git
   ```

2. Installation:

   Nf-NeST comes with a Doker image that can be run pipeline with virtual environment. To setup up Docker image, run the following command from the Nf-NeST directory. 

   ```
   cd  Nf-NeST
   docker pull supark87/nfnest:latest
   ``` 
3. Input raw sequences name
   
   Paired reads should be separated including 'R1' and 'R2' in sequence name. 

<a id="First"></a>
## Your first analysis

   NeST was conceptualized to identify mutations that confer anti-malarial drug resistance in *P.falciparum* (Talundzic et al., 2018). It was also applied for the detection of antibiotic drug resistance in *M.tubercolosis* (Colman et al., 2015). 

1. Test Run
   - Use Docker
   10 samples of *P.falciparum* from NCBI are located in /testrun/fastq folder.
   To execute Nf-NeST pipeline on this samples, run this command line
   ```
   docker run -v $(pwd)/testrun:/data/testrun -v $(pwd)/pyscripts:/data/pyscripts -ti supark87/nfnest:latest ./nextflow run nfNeST_ver02.nf -c   ./testrun/nextflow1.config -with-report ./testrun/test_output.html
   
   ```
- Use singularity (scicomp)
   ```
   module load singularity/3.5.3
   singularity pull docker://supark87/nfnest_singularity
   singularity run -B $(pwd)/testrun:/data/testrun -B $(pwd)/pyscripts:/data/pyscripts nfnest_singularity_latest.sif nextflow run nfNeST_ver02.nf -c /data/testrun/nextflow1.config -with-report /data/testrun/test_output.html

   ```


   ```
## Executing your own analysis using Nf-NeST:   
      
   Copy your inputs under the folder /inputfiles/. By default, configuration file for this folder is in here as 'nextflow.config'
      
   Nf-NeST can be executed on your own dataset using the following command:

      ```
      docker run -v $(pwd)/inputfiles:/data/inputfiles -v $(pwd)/pyscripts:/data/pyscripts-ti supark87/nfnest:latest ./nextflow run nfNeST_ver02.nf \
      -c ./inputfiles/nextflow.config -with-report ./inputfiles/output/output.html
      ```

      The details about the required input formats are listed in the next section.

<a id="inputs"></a>
## Input standardization:

NeST is designed to reduce the amount of user intervention with regards to inputs that the user needs to provide. However to enable standardization of inputs across all organisms we require that a particular file format be followed for the three inputs listed below:

1. Fastq files:

   The PrepInputs module in NeST highly simplifies the management of fastq files. The module accepts two input formats.
   - Input directory path:

      This just requires the user to provide the path to a folder containing fastq files. The files are recognized by the file extension, so the files must have either ```fq```, ```fq.gz```, ```fastq``` or ```fastq.gz``` file extensions. The name convention of paired file can be ```_R1```.


2. BED format:

   The BED (Browser Extensible Data) is an easy and lightweight format to list annotations for a genome. NeST uses a full BED or BED 12 column format file as a guide to annotate variants with codon and amino acid changes. The example file listed below shows the details of how to structure the BED file. The separation of contig, gene and exon level information makes this format highly portable across genomes. The BED 12 column format for most organisms can be export from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). A detail explanation of the BED format can be found [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

   ```
   #contig start stop gene score strand  CDSstart  CDSstop rbg NoOfExons ExonLengths ExonStarts
   PfCRT 1	3399	PfCRT	.	+	95	3191	0	13	90,268,172,132,71,75,82,50,56,92,44,54,76,	96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115,
   MT	1	5967	COXIII	.	-	734	1573	0	1	839,	734,
   MT	1	5967	COL	.	+	1933	3471	0	1	1538,	1933,
   MT	1	5967	CYTOb	.	+	3492	4622	0	1	1130,	3492,
   PfDHFR	1	1827	PfDHFR	.	+	1	1827	0	1	1827,	1,
   PfDHPS	1	2417	PfDHPS	.	+	1	2417	0	3	135,1868,115,	1,313,2302,
   PfK13	1	2181	PfK13	.	+	1	2181	0	1	2181,	1,
   PfMDR1	1	4260	PfMDR1	.	+	1	4260	0	1	4260,	1,
   ```

3. Variant of Interest:

   The Summarize module within NeST, allows for easy summarization of variants called from all samples in a study. If a user specifies a list of variants of interest, a separate table will be created for these set of variants. The variants can be specified in ```.tsv```, ```.csv```, ```.xlsx``` format. And follows the format listed below

   | Chrom  | Gene   | RefAA | AAPos | AltAA |
   |:------:|:------:|:-----:|:-----:|:-----:|
   | PfCRT  | PfCRT  |   C   |   72  |   S   |
   | PfCRT  | PfCRT  |   V   |   73  |   V   |
   | PfMDR1 | PfMDR1 |   N   |   86  |   Y   |
   | PfMDR1 | PfMDR1 |   Y   |   184 |   F   |
   | MT     | CYTOb  |   I   |   258 |   M   |



<a id="outputs"></a>
## Output Description

Output file folder will be created under /inputfiles/

Output table will be found in /inputfiles/output/visualization/ 

SNP information for each sample with separate file could be found in /inputfiles/output/snpfilter

1. Report files:

   Nf-NeST produces table reports that summarize the different types of variants found in the sample. All the tables will be stored under the output directory. The table below describes the different files that are generated by Nf-NeST.

   |                   File                    |                                          Description                                                                                                      |
   |:------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------|
   | summary.csv                  | This file contains the calls for each of the variants of interst, for each of the samples. It provides, allele frequency, mutation type, and amino acid sequence at the position of interest |
   | haplotype.csv | This file lists the haplotype(pseudo-haplotype) that were constructed from alleles that have allele frequency greater or equal to 0.95         |
   | PowerBI_input.csv            | This file will be fed into PowerBI tool as inputfile for interactive visualization                |
   | files under output/combined_json                       | JSON file with sample meta information and variant calls for the all samples in the study                                                                 |
   | output_report*.html          | This file would be created in output directory. You can see whole monitoring process there with failed jobs with time, cpu, and failed reports as needed.                                    |

2. PowerBI reports 

   Nf-NeST produces PowerBI input files that would be fed into PowerBI interactive visualization tool.

   - Under development

    

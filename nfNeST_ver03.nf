#!/usr/bin/env nextflow
params.reads = params.input.fastq_path+'*{R1,R2}*.fastq.gz'
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
bbduk_path = "$params.input.bbduk_def"
params.genome = "$baseDir/ref/pfalciparum/New_6_genes.fa"
params.bed = "$baseDir/ref⁩/pfalciparum⁩/New_6_genes.bed"
params.fas = "$baseDir/ref/pfalciparum/adapters.fa"
params.gatk= "$baseDir/gatk-4.1.9.0/gatk"
mode = "$params.input.mode"
vcfmode = "$params.input.vcfmode"
logpath="$baseDir/.nextflow.log"
workpath="$baseDir/work"
CTexcelpath="$baseDir/ANG_2019_TES_master.xlsx"
datasettype = "$params.input.datasettype"
configfiles="$baseDir/6genes_ver3"
pyscripts="$baseDir/pyscripts"
//snpfilter_scripts="$baseDir/pyscript"




process combineFastq {
    publishDir "$params.output.folder/trimFastq/${pair_id}", pattern: "*.fastq", mode : "copy"
    input:
        set val(pair_id), path(fastq_group) from fastq_path
    output:
        tuple val(pair_id), path("${pair_id}_R1.fastq.gz"), path("${pair_id}_R2.fastq.gz") into comb_out

    script:
        """
        cat *R1*.fastq.gz > ${pair_id}_R1.fastq.gz
        cat *R2*.fastq.gz > ${pair_id}_R2.fastq.gz
        """
}

comb_out.into{preqc_path; trim_path}

process preFastQC {
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process trimFastq {
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    publishDir "$params.output.folder/trimFastq/${sample}/Stats", pattern: "*.txt", mode : "copy"
    errorStrategy 'ignore'

    input:
        set val(sample), path(read_one), path(read_two) from trim_path
        path fas from params.fas
    
    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """          
        $baseDir/bbmap/bbduk.sh -Xmx1g ktrimright=t k=27 hdist=1 edist=0 ref=${fas} \\
        qtrim=rl trimq=35 minlength=150 trimbyoverlap=t minoverlap=24 ordered=t qin=33 in=${read_one} in2=${read_two} \\
        out=${sample}_trimmed_R1.fastq out2=${sample}_trimmed_R2.fastq stats=${sample}_stats.txt 
        """
}

trim_out.into{align_path; postqc_path}

process postFastQC {
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from postqc_path
    
    output:
        tuple val(sample), path("*") into postqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process buildIndex {
    tag "$genome.baseName"
    publishDir "$params.output.folder/Bowtie2Index", mode : "copy"

    input:
    path genome from params.genome

    output:
    file 'genome.index*' into index_ch

    script:
    if( mode == 'Bowtie')
    """
    bowtie2-build ${genome} genome.index
    """
}



process alignReads {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from align_path
        file index from index_ch
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out2

        tuple val(sample), path("${sample}.sam") into sam_out

    script:
        index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
        if( mode == 'Bowtie')
            """
            bowtie2 --very-sensitive  --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
            -x $baseDir/$params.output.folder/Bowtie2Index/$index_base \\
            -1 ${read_one} -2 ${read_two} -p 4 -S ${sample}.sam

            """
        else if( mode == 'Bwa' )
            """
            bwa mem -t 4 ${genome} ${read_one} ${read_two} > ${sample}.sam 
            """
        else if( mode == 'BBMap' )
            """
            bbmap ref=${genome} in=${read_one} in2=${read_two} out=${sample}.sam 
            """
        else if( mode == 'Snap' )
            """
            snap paired ${genome} ${read_one} ${read_two} -t 4 -o -sam ${sample}.sam $task.cpus \\
            """
}

process processAlignments {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.ba*", mode : "copy"
    errorStrategy 'retry'

    input:
        set val(sample), path(sam_path), path(align_met) from align_out
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}_SR.bam") into postal_out
        tuple val(sample), path("${sample}_SR.bam") into postal_out2
        tuple val(sample), path("${sample}_SR.bam") into postal_out3
        tuple val(sample), path("${sample}_SR.bam") into postal_out4
        tuple val(sample), path("${sample}_SR.bam") into postal_out5
        tuple val(sample), path("${sample}_SR.bam") into postal_out6


    script:
        """
        java -jar $baseDir/picard.jar AddOrReplaceReadGroups -I ${sam_path} -O ${sample}_SR.bam -SORT_ORDER coordinate --CREATE_INDEX true \\
        -LB ExomeSeq -DS ExomeSeq -PL Illumina -CN AtlantaGenomeCenter -DT 2016-08-24 -PI null -ID ${sample} \\
        -PG ${sample} -PM ${sample} -SM ${sample} -PU HiSeq2500 > file
        """
    
}

process getcoverage {

    container 'zlskidmore/samtools'
    publishDir "$params.output.folder/samtoolcoverage/${sample}", mode : "copy"
    
    input:
    set val(sample), path(bam_path) from postal_out6 
    tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") from align_out2

    output:
    file "${sample}_coverage.txt" into coverage
    tuple val(sample), path("${sample}_finalcoverageall.txt") into finalcoverage
    tuple val(sample), path("${sample}_finalcoverageall.txt") into finalcoverage2
    tuple val(sample), path("${sample}_depth.txt") into depth
    script:
    
    """
   $baseDir/samtools-1.13/samtools index ${bam_path} 
   $baseDir/samtools-1.13/samtools coverage  ${bam_path} -o ${sample}_coverage.txt    
   $baseDir/samtools-1.13/samtools depth  ${bam_path} -a -o ${sample}_depth.txt

   cat *coverage.txt > ${sample}_finalcoverageall.txt
   cat *all.txt > coverage_combine.txt
   cat *depth.txt > depth_combine.txt

    """
}
process GenerateVCF {
    publishDir "$params.output.folder/GenerateVCFSam/${sample}", pattern: "*.vcf", mode : "copy"
    input:
        set val(sample), path(bam_path) from postal_out
        path genome from params.genome

    output:
        tuple val(sample), path("${sample}-1.vcf"), path("${sample}-2.vcf"), path("${sample}-3.vcf") into vcf_out1

    script:
        """
        $baseDir/bcftools-1.9/bcftools mpileup -f ${genome} ${bam_path} > ${sample}.mpileup
        $baseDir/bcftools-1.9/bcftools call -vm ${sample}.mpileup  > ${sample}-1.vcf
        $baseDir/gatk-4.2.0.0/gatk HaplotypeCaller -R $baseDir/ref/pfalciparum/New_6_genes.fa -I ${bam_path} -O ${sample}-2.vcf --base-quality-score-threshold 6 --max-reads-per-alignment-start 0 --min-base-quality-score 0 --standard-min-confidence-threshold-for-calling 0 
        freebayes -f ${genome} -F 0.01 -E 3 --report-all-haplotype-alleles --haplotype-length 1  ${bam_path} > ${sample}-3.vcf

        """
                //freebayes -f ${genome} -F 0.01 -E 3 -t $baseDir/6genes_ver3/6genes_ver3/targe_.bed --report-all-haplotype-alleles ${bam_path} > ${sample}-3.vcf

}


process annotate {
    publishDir "$params.output.folder/annotate/${sample}", mode : "copy"
    errorStrategy 'ignore'

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(vcf_path3) from vcf_out1
        path(configpath) from configfiles


    output:
        tuple val(sample), path("${sample}-1.ann.vcf"), path("${sample}-2.ann.vcf"),path("${sample}-3.ann.vcf") into anno_out1 

    script:
        """
       java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar build -c ${configpath}/snpEff.config -v 6genes_ver3
       java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config 6genes_ver3 -interval $baseDir/ref/pfalciparum/annotations2.bed ${vcf_path1} > ${sample}-1.ann.vcf
       java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config 6genes_ver3 -interval $baseDir/ref/pfalciparum/annotations2.bed ${vcf_path2} > ${sample}-2.ann.vcf
       java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config 6genes_ver3 -interval $baseDir/ref/pfalciparum/annotations2.bed ${vcf_path3} > ${sample}-3.ann.vcf

        """
        //       java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar build -c ${configpath}/snpEff.config -v 6genes_ver3

        // java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config New_6_genes -interval $baseDir/ref/pfalciparum/annotation.bed ${vcf_path1} > ${sample}-1.ann.vcf
        // java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config New_6_genes  -interval $baseDir/ref/pfalciparum/annotation.bed ${vcf_path2} > ${sample}-2.ann.vcf
        // java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config New_6_genes -interval $baseDir/ref/pfalciparum/annotation.bed ${vcf_path3} > ${sample}-3.ann.vcf

}

process vartpype {
    publishDir "$params.output.folder/vartype/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(vcf_path3) from anno_out1

    output:
        tuple val(sample), path("${sample}_1.vartype.vcf"), path("${sample}_2.vartype.vcf"), path("${sample}_3.vartype.vcf") into vartype_out1

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar varType ${vcf_path1} > ${sample}_1.vartype.vcf
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar varType ${vcf_path2} > ${sample}_2.vartype.vcf
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar varType ${vcf_path3} > ${sample}_3.vartype.vcf
        """
}

process merge {
    publishDir "$params.output.folder/final_vcf/${sample}", mode : "copy"
    errorStrategy 'ignore'

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(vcf_path3),path(bam_path) from vartype_out1.join(postal_out4)
        path(pyscripts_path) from pyscripts

    output:
        tuple val(sample), path("final_${sample}.vcf") into merge_out
        tuple val(sample), path("final_${sample}.vcf") into merge_out2
        tuple val(sample), path("final_${sample}.vcf") into merge_out3


    script:
        """
        samtools index ${bam_path}
        python3 ${pyscripts_path}/annotate2.py -r $baseDir/ref/pfalciparum/New_6_genes.fa -b $baseDir/ref/pfalciparum/annotation2.bed -o ${sample} -v1 ${vcf_path1} -v2 ${vcf_path2} -v3 ${vcf_path3} -m ${bam_path} -voi $baseDir/ref/pfalciparum/voinew2.csv -name ${sample}
        """
}

process filter {
    publishDir "$params.output.folder/filter/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from merge_out

    output:
        tuple val(sample), path("${sample}_filtered.vcf") into filter_out

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar filter -f ${vcf_path1} " ( VARTYPE = 'SNP'  ) " > ${sample}_filtered.vcf
        """
}

process extract {
    publishDir "$params.output.folder/extract/${sample}", mode : "copy"
    errorStrategy 'ignore'
    
    input:
        set val(sample), path(vcf_path1) from filter_out

    output:
        tuple val(sample), path("final_${vcf_path1}ext.vcf") into extract_out

    script:
    """
    java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar extractFields ${vcf_path1} CHROM POS REF ALT VARTYPE Confidence Sources DP4 DP AD GEN[*].AD "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > final_${vcf_path1}ext.vcf
    """
}

process spread {
    publishDir "$params.output.folder/spread/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path) from extract_out
        path(pyscripts_path) from pyscripts

    output:
        tuple val(sample), path("fixedPOS${vcf_path}") into spread_out

    script:
    """
    python3 ${pyscripts_path}/vcfcsv3.py -n ${vcf_path}
    """
}

process snpfilter {
    publishDir "$params.output.folder/snpfilter/${sample}", pattern: "*.csv", mode : "copy"
    publishDir "$params.output.folder/log/", pattern: "*.txt", mode : "copy"
    publishDir "$params.output.folder/forjson/", pattern: "*.csv", mode : "copy"

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(bam_path1) from spread_out.join(merge_out2).join(postal_out5)
        path(pyscripts_path) from pyscripts

    output:
        tuple val(sample), path("${sample}.csv") into snpfilter_out
        tuple val(sample), path("${sample}.csv") into snpfilter_out2


    script:
        """
        python3 ${pyscripts_path}/main.py -v1 ${vcf_path2} -v2 ${vcf_path1} -b1 ${bam_path1} -o1 ${sample} -e1 $baseDir/ref/pfalciparum/candidates.xlsx -e2 $baseDir/ref/pfalciparum/voinew4.csv -f1 $baseDir/ref/pfalciparum/New_6_genes.fa -b2 $baseDir/ref/pfalciparum/mdr.bed
        """
}
        //python3 ${pyscripts_path}/snpreport1-54-VAF_SP.py -v1 ${vcf_path2} -v2 ${vcf_path1} -b1 ${bam_path1} -b2 $baseDir/ref/pfalciparum/mdr.bed -o1 ${sample} -e1 $baseDir/ref/pfalciparum/candidates.xlsx -e2 $baseDir/ref/pfalciparum/voinew4.csv -f1 $baseDir/ref/pfalciparum/New_6_genes.fa

        //mv $logfile "nextflowlog.txt"



process tojson {
    publishDir "$params.output.folder/combinedjson/", pattern: "*.json", mode : "copy"
    input:
        set val(sample), path("${sample}.csv") from snpfilter_out
        path(pyscripts_path) from pyscripts

    output:
        tuple val(sample), path("wholecombine.json") into tojson_out
    script:
        """
        python3 ${pyscripts_path}/tojson2.py -d1 $baseDir/$params.output.folder/forjson
        """
}

process getcutoff{

    publishDir "$params.output.folder/cutoff/", mode : "copy"
    errorStrategy 'ignore'


    input:
    file("*") from finalcoverage.collect()
    file("*") from depth.collect()
    path(pyscripts_path) from pyscripts

    output:
    file("first_q_coverage.csv")
    file("first_q_depth.csv")     
    file("list.csv") into list
    file("list1.csv") into list1

    script:
    """
    python3 $pyscripts/getcutoff.py -c $baseDir/$params.output.folder/samtoolcoverage
    """
 
}



process viz {   
    publishDir "$params.output.folder/visualization/", mode : "copy"
    errorStrategy 'ignore'

    input:
        path(logfile) from logpath
        path(pyscripts_path) from pyscripts
        file("*") from snpfilter_out2.collect()
        file(list) from list
        file(list1) from list1
      
    output:
        file("summary_reportable.csv")
        file("NCBI_feature_column.csv")
        file("PowerBI_input.csv")
        file("cutoff_SNP.csv")
        file("summary.csv")
        file("haplotype.csv")
        file("nextflowlog.txt") into logpath1
        file("haplotype_count.csv")

    script:
        """
        python3 ${pyscripts_path}/visualization_final2.py -voi $baseDir/ref/pfalciparum/voinew4.csv -cod $baseDir/$params.output.folder/cutoff/
        mv $logfile  "nextflowlog.txt"
        """
} 

process errortrack {
    publishDir "$params.output.folder/errortrack/", mode : "copy"
    errorStrategy 'ignore'
    input:

        file(logfile) from logpath1

        path(work) from workpath
        path(pyscripts_path) from pyscripts

        path(CT) from CTexcelpath
    output:
        file("out.csv")
    script:
        if( datasettype == 'Angola')
            """
            python3 ${pyscripts_path}/errortrackn6.py -n1 ${logfile} -w1 ${work} -x1 ${CT}
            """
        else
            """
           python3 ${pyscripts_path}/errortrackn7.py -n1 ${logfile} -w1 ${work}
            """
}

// process sendemail{
//     """
//     sendMail {
//     to 'nej1@cdc.gov'
//     from 'nej1@cdc.gov'
//     attach "$params.output.folder/visualization/PowerBI_input.csv"
//     attach "$params.output.foler/snpfilter/*/*.csv"
//     subject 'Your nfNeST job is done'

//     '''
//     Hi there,
//     Please find the files attached for nfNeST results
//     Thanks,
//     Subin
//     '''
    
//     }
//     """
// }

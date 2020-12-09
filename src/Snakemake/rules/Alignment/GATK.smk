
chrom_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

rule GATK_recal_step1:
    input:
        bam = "bam/{sample}-sort-cumi.bam",
        bai = "bam/{sample}-sort-cumi.bam.bai",
        dbsnp = "/data/ref_genomes/hg19/variation/dbsnp_138.vcf.gz", #config
        bed = config["bed"]["bedfile"],
        ref = config["reference"]["ref"]
    output:
        grp = "bam/{sample}-sort-cumi-recal.grp"
    params:
        "--interval_set_rule INTERSECTION -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment"
    log:
        "logs/gatk3/recal_step1_{sample}.log"
    singularity:
        config["singularity"]["gatk3"]
    threads:
        10
    shell:
        "(java -jar -Xms1000m -Xmx50960m /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -nct {threads} -I {input.bam} -o {output.grp} -R {input.ref} --knownSites {input.dbsnp} -L {input.bed} {params}) &> {log}"


rule GATK_recal_step2:
    input:
        grp = "bam/{sample}-sort-cumi-recal.grp",
        bam = "bam/{sample}-sort-cumi.bam",
        bai = "bam/{sample}-sort-cumi.bam.bai",
        ref = config["reference"]["ref"]
    output:
        bam = "bam/{sample}-sort-cumi-recal.bam"
    params:
        "-jdk_deflater -jdk_inflater -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment"
    log:
        "logs/gatk3/recal_step2_{sample}.log"
    singularity:
        config["singularity"]["gatk3"]
    threads:
        2
    shell:
        "(java -jar -Xms1000m -Xmx91728m /usr/GenomeAnalysisTK.jar -T PrintReads -nct {threads} -R {input.ref} -I {input.bam} -BQSR {input.grp} -o {output.bam} {params}) &> {log}"


rule Split_bam_realign:
    input:
        bam = "bam/{sample}-sort-cumi-recal.bam",
        #bai = "DNA_bam/{sample}-ready.bam.bai"
        # vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    output:
        bam = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam",
        bai = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam.bai"
    log:
        "logs/gatk3/split_bam_realign_{sample}-sort-cumi-recal-{chr}.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools view -b {input.bam} {wildcards.chr} > {output.bam} && samtools index {output.bam}) &> {log}"


rule GATK_realign_step1:
    input:
        bam = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam",
        bai = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam.bai",
        ref = config["reference"]["ref"],
        indels = "/data/ref_genomes/hg19/variation/Mills_and_1000G_gold_standard.indels.vcf.gz" #config
    output:
        intervals = "bam/{sample}-sort-cumi-recal-realign.{chr}.intervals"
    params:
        "--interval_set_rule INTERSECTION -L {chr} -l INFO -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment"
    log:
        "logs/gatk3/realign_step1_{sample}_{chr}.log"
    singularity:
        config["singularity"]["gatk3"]
    shell:
        "(java -jar -Xms500m -Xmx3500m /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {input.ref} -I {input.bam} --known {input.indels} -o {output.intervals} {params}) &> {log}"


rule GATK_realign_step2:
    input:
        bam = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam",
        bai = "bam/realign_temp/{sample}-sort-cumi-recal.{chr}.bam.bai",
        ref = config["reference"]["ref"],
        indels = "/data/ref_genomes/hg19/variation/Mills_and_1000G_gold_standard.indels.vcf.gz", #config
        intervals = "bam/{sample}-sort-cumi-recal-realign.{chr}.intervals",
    output:
        bam = "bam/realign_temp/{sample}-sort-cumi-recal-realign.{chr}.bam"
    params:
        "-L {chr} -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment"
    log:
        "logs/gatk3/realign_step2_{sample}_{chr}.log"
    singularity:
        config["singularity"]["gatk3"]
    shell:
        "(java -jar -Xms909m -Xmx6363m /usr/GenomeAnalysisTK.jar -T IndelRealigner -R {input.ref} -I {input.bam} --targetIntervals {input.intervals} --knownAlleles {input.indels} -o {output.bam} {params}) &> {log}"


rule Merge_bam_gatk3:
    input:
        bams = expand("bam/realign_temp/{{sample}}-sort-cumi-recal-realign.{chr}.bam", chr=chrom_list)
    output:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai"
    log:
        "logs/gatk3/merge_bam_{sample}.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools merge {output.bam} {input.bams} && samtools index {output.bam}) &> {log}"

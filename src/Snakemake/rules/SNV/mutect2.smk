

localrules: split_bedfile, fixSB
chrom_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

rule Split_bam:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai"
        # vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    output:
        bam = temp("mutect2/bam_temp/{sample}-ready.{chr}.bam"),
        bai = temp("mutect2/bam_temp/{sample}-ready.{chr}.bam.bai")
    log:
        "logs/split_bam_{sample}-ready-{chr}.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools view -b {input.bam} {wildcards.chr} > {output.bam} && samtools index {output.bam}) &> {log}"

rule split_bedfile:
    input:
        config["bed"]["bedfile"]
    output:
        temp("mutect2/bedfile.{chr}.bed")
    log:
        "logs/variantCalling/split_bed.{chr}.log"
    shell:
        "(grep -w {wildcards.chr} {input} > {output}) &> {log}"


rule Mutect2_gatk3:
    input:
        bam = "mutect2/bam_temp/{sample}-ready.{chr}.bam",
        bai = "mutect2/bam_temp/{sample}-ready.{chr}.bam.bai",
        fasta = config["reference"]["ref"],
        bed = "mutect2/bedfile.{chr}.bed"
    output:
        bam = temp("mutect2/bam_temp2/{sample}-ready.{chr}.indel.bam"),
        bai = temp("mutect2/bam_temp2/{sample}-ready.{chr}.indel.bai"),
        #stats = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf.gz.stats"),
        #vcf = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf.gz")
        vcf = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf")
        #vcf_tbi = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf.gz.tbi")
    params:
        annotation = "--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation GCContent --annotation HaplotypeScore --annotation HomopolymerRun --annotation DepthPerAlleleBySample --annotation Coverage",
        parameters = "--interval_set_rule INTERSECTION -drf DuplicateRead -ploidy 2 -U LENIENT_VCF_PROCESSING ",
        filter = "--read_filter BadCigar --read_filter NotPrimaryAlignment"
    log:
        "logs/variantCalling/mutect2_{sample}.{chr}.log"
    singularity:
        config["singularity"]["mutect2_gatk3"]
    shell:
        "(java -jar -Xmx4g /usr/GenomeAnalysisTK.jar -T MuTect2 {params.annotation} {params.parameters} {params.filter} -R {input.fasta} -I:tumor {input.bam} -L {input.bed} --bamOutput {output.bam} -o {output.vcf}) &> {log}"


# rule Mutect2:
#     input:
#         bam = "mutect2/bam_temp/{sample}-ready.{chr}.bam",
#         bai = "mutect2/bam_temp/{sample}-ready.{chr}.bam.bai",
#         fasta = config["reference"]["ref"],
#         bed = "mutect2/bedfile.{chr}.bed"
#     output:
#         bam = temp("mutect2/bam_temp2/{sample}-ready.{chr}.indel.bam"),
#         bai = temp("mutect2/bam_temp2/{sample}-ready.{chr}.indel.bai"),
#         stats = temp("mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz.stats"),
#         vcf = temp("mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz"),
#         vcf_tbi = temp("mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz.tbi")
#     log:
#         "logs/variantCalling/mutect2_{sample}.{chr}.log"
#     singularity:
#         config["singularity"]["mutect2"]
#     shell:
#         "(gatk --java-options '-Xmx4g' Mutect2 -R {input.fasta} -I {input.bam} -L {input.bed} --bam-output {output.bam} -O {output.vcf}) &> {log}"
#
# rule filterMutect2:
#     input:
#         vcf = "mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz",
#         vcf_tbi = "mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz.tbi",
#         stats = "mutect2/bam_temp2/{sample}.{chr}.mutect2.unfilt.vcf.gz.stats",
#         fasta = config["reference"]["ref"]
#     output:
#         vcf = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf.gz"),
#         vcf_tbi = temp("mutect2/filteringStats/{sample}.{chr}.mutect2.vcf.gz.tbi")
#     log:
#         "logs/variantCalling/mutect2/filter_{sample}.{chr}.log"
#     singularity:
#         config["singularity"]["mutect2"]
#     shell:
#         "(gatk --java-options '-Xmx4g' FilterMutectCalls --max-alt-allele-count 3 --max-events-in-region 8 -R {input.fasta} -V {input.vcf} -O {output.vcf}) &> {log}"

rule Merge_vcf:
    input:
        #vcf = expand("mutect2/filteringStats/{{sample}}.{chr}.mutect2.vcf.gz", chr=chrom_list)
        vcf = expand("mutect2/filteringStats/{{sample}}.{chr}.mutect2.vcf", chr=chrom_list)
    output:
        temp("mutect2/{sample}.mutect2.SB.vcf")
    log:
        "logs/variantCalling/mutect2/merge_vcf_{sample}.log"
    singularity:
        config["singularity"]["bcftools"]
    shell:
        "(bcftools concat -o {output} -O v {input} ) &> {log}"

rule fixSB:
    input:
        "mutect2/{sample}.mutect2.SB.vcf"
    output:
        temp(touch("mutect2/{sample}.SB.done"))
    log:
        "logs/variantCalling/mutect2/{sample}.fixSB.log"
    shell:
        "(sed -i 's/=SB/=SB_mutect2/g' {input}  && sed -i 's/:SB/:SB_mutect2/g' {input}) &> {log}"

rule mutect2HardFilter:
    input:
        vcf = "mutect2/{sample}.mutect2.SB.vcf",
        wait = "mutect2/{sample}.SB.done"
    output:
        temp("mutect2/{sample}.mutect2.fixAF.vcf")
    log:
        "logs/variantCalling/mutect2/{sample}.hardFilt.log"
    singularity:
        config["singularity"]["python"]
    shell:
        "(python3.6 src/hardFilter_PASS_mutect2.py {input.vcf} {output}) &> {log}"

rule Merge_bam:
    input:
        bams = expand("mutect2/bam_temp2/{{sample}}-ready.{chr}.indel.bam", chr=chrom_list)#["mutect2/perChr/{sample}" + str(c) + ".indel.bam" for c in chrom_list]
    output:
        bam = "DNA_bam/mutect2_bam/{sample}-ready.indel.bam",
        bai = "DNA_bam/mutect2_bam/{sample}-ready.indel.bam.bai"
    log:
        "logs/variantCalling/merge_bam_{sample}.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools merge {output.bam} {input.bams} && samtools index {output.bam}) &> {log}"


rule vardict:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        ref = config["reference"]["ref"],
        bed = config["bed"]["bedfile"]
    output:
        temp("vardict/{sample}.vardict.vcf")
    params:
        vardict_singularity = config["singularity"]["execute"] + config["singularity"]["vardict"],
        bcftools_singularity = config["singularity"]["execute"] + config["singularity"]["bcftools"],
        af = "0.01"
    log:
        "logs/variantCalling/vardict/{sample}.log"
    threads:
        5
    shell:
        "({params.vardict_singularity} /opt/VarDict-1.5.1/bin/VarDict -G {input.ref} -f {params.af} -I 200 -th {threads} -N '{wildcards.sample}' -z -c 1 -S 2 -E 3 -g 4 -Q 10 -F 0x700 -b {input.bam} {input.bed} |"
        " {params.vardict_singularity} /opt/VarDict/teststrandbias.R |"
        " {params.vardict_singularity} /opt/VarDict/var2vcf_valid.pl -A -N '{wildcards.sample}' -E -f {params.af} |"
        " {params.bcftools_singularity} bcftools filter -i 'QUAL >= 0' |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $5) }} {{print}}' |"
        " awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {{next}} {{print}}' > {output}) 2> {log}"


rule sortVardict:
    input:
        "vardict/{sample}.vardict.vcf"
    output:
        "vardict/{sample}.vardict.okAF.vcf"
    singularity:
        config["singularity"]["bcftools"]
    log:
        "logs/variantCalling/vardict/{sample}.sort.log"
    shell:
        "(bgzip {input} && tabix {input}.gz && "
        "bcftools sort -o {output} -O v {input}.gz) &> {log}"

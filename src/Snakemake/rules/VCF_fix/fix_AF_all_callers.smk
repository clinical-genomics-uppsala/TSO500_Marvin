
rule fixAF:
    input:
        "{method}/{sample}.{method}.fixAF.vcf"
    output:
        temp("{method}/{sample}.{method}.okAF.vcf")
    log:
        "logs/variantCalling/fixAF/{method}/{sample}.log"
    singularity:
        config["singularity"]["python"]
    shell:
        "(python3.6 src/fix_af.py {input} {output}) &> {log}"

localrules: bgzipCallers
rule bgzipCallers:
    input:
        vcf = "{method}/{sample}.{method}.okAF.vcf" #[m+"/{sample}."+m+".vcf" for m in config["methods"]]
    output:
        #vcf = temp("{method}/{sample}.{method}.vcf.gz"),
        #tabix = temp("{method}/{sample}_{seqID}.{method}.vcf.gz.tbi")
        vcf = "{method}/{sample}.{method}.okAF.vcf.gz",
        tabix = "{method}/{sample}.{method}.okAF.vcf.gz.tbi"
    log:
        "logs/variantCalling/bgzip/{method}/{sample}.log"
    singularity:
        config["singularity"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) 2> {log}"

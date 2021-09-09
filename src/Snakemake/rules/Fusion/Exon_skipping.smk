
rule Exon_skipping:
    input:
        bed=config["bed"]["MET"],
        junction="STAR/{sample}SJ.out.tab",
    output:
        results="Results/RNA/{sample}/Fusions/{sample}_Met_exon_skipping.txt",
    log:
        "logs/Fusion/{sample}_exon_skipping.log"
    singularity:
        config["singularity"]["python"]
    shell:
        "(python3.6 src/Exon_skipping.py {input.bed} {input.junction} {output.results}) &> {log}"

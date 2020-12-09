
rule fastqc_bam:
    input:
        #"fastq/{sample}_R1-fastq.gz" ##one for each R1 and one for R2 should be from a samples.yaml file
        #"DNA_bam/{sample}-ready.bam"
        "bam/{sample}-sort.bam"
    output:
        html="qc/{sample}/{sample}-sort_fastqc.html",
        zip="qc/{sample}/{sample}-sort_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        outdir = "qc/{sample}/"
    log:
        "logs/qc/fastqc/{sample}-sort.log"
    threads: 10
    singularity:
        config["singularity"]["fastqc"]
    shell:
        "(fastqc --quiet -t {threads} --outdir {params.outdir} {input}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/fastqc"

# rule fastqcR2:
#     input:
#         "data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq" ##one for each R1 and one for R2 should be from a samples.yaml file
#     output:
#         html="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.html",
#         zip="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params: ""
#     log:
#         "logs/qc/fastqc/{sample}_{seqID}_R2_trimmed.log"
#     singularity:
#         config["singularitys"]["fastqc"]
#     wrapper:
#         "0.38.0/bio/fastqc"

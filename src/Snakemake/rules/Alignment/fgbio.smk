
rule fgbio:
    input:
        bam = "bam/{sample}-sort.bam",
        ref = config["reference"]["ref"]
    output:
        fq1 = temp("fastq_temp/{sample}-cumi-R1.fq.gz"),
        fq2 = temp("fastq_temp/{sample}-cumi-R2.fq.gz"),
        qc = "qc/{sample}/{sample}_fgbio.txt"
    params:
        bam_tmp = "bam/{sample}-cumi-1-bamtofastq-tmp",
        fgbio_singularity = config["singularity"]["execute"] + config["singularity"]["fgbio"],
        bamtofastq_singularity = config["singularity"]["execute"] + config["singularity"]["bamtofastq"]
    log:
        "logs/fgbio/{sample}.log"
    shell:
        "{params.fgbio_singularity} fgbio GroupReadsByUmi -i {input.bam} --edits=1 --min-map-q=1 -t RX -s adjacency -f {output.qc}"
        " | {params.fgbio_singularity} fgbio CallMolecularConsensusReads -i /dev/stdin -o /dev/stdout"
        " --min-input-base-quality=2 --min-reads=1 --max-reads=100000 --output-per-base-tags=false --sort-order=:none:"
        " | {params.fgbio_singularity} fgbio FilterConsensusReads -i /dev/stdin -o /dev/stdout -r {input.ref} --min-reads=1 --min-base-quality=13"
        " | {params.bamtofastq_singularity} bamtofastq collate=1 tags=cD,cM,cE gz=1 T={params.bam_tmp} F={output.fq1} F2={output.fq2}"
        " > {log}"

rule bwa_mem_fgbio:
    input:
        reads = ["fastq_temp/{sample}-cumi-R1.fq.gz", "fastq_temp/{sample}-cumi-R2.fq.gz"]
    output:
        #"DNA_bam/{sample}-ready.bam"
        "bam/{sample}-sort-cumi.bam"
    log:
        "logs/map/bwa/{sample}.log"
    params:   #-M
        index = config["reference"]["ref"],
        extra = r"-C -c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
        # sort = "samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order = "coordinate"  # Can be 'queryname' or 'coordinate'.
        # sort_extra = ""            # Extra args for samtools/picard.
    threads: 10
    singularity:
        config["singularity"]["bwa"]
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -@ {threads} -m 3G -o {output} - ) &> {log}"


rule samtools_index_fgbio:
    input:
        #"DNA_bam/{sample}-ready.bam"
        "bam/{sample}-sort-cumi.bam"
    output:
        #"DNA_bam/{sample}-ready.bam.bai"
        "bam/{sample}-sort-cumi.bam.bai"
    log:
        "logs/map/samtools_index/{sample}-ready.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools index {input} {output}) &> {log}"


rule bwa_mem:
    input:
        reads = ["fastq/DNA/{sample}_R1.fastq.gz", "fastq/DNA/{sample}_R2.fastq.gz"]
    output:
        "bam/{sample}-sort.bam"
    log:
        "logs/map/bwa/{sample}.log"
    params:   #-M
        bwa_singularity = config["singularity"]["execute"] + config["singularity"]["bwa"],
        bamsormadup_singularity = config["singularity"]["execute"] + config["singularity"]["bamsormadup"],
        umis_singularity = config["singularity"]["execute"] + config["singularity"]["umis"],
        samtools_singularity = config["singularity"]["execute"] + config["singularity"]["samtools"],
        index = config["reference"]["ref"],
        extra = r"-c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
        # sort = "samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order = "coordinate", # Can be 'queryname' or 'coordinate'.
        # sort_extra = ""            # Extra args for samtools/picard.
        tmp_dir = "tmpfile=bam/{sample}"
    threads: 10
    #singularity:
    #    "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "({params.bwa_singularity} bwa mem -t {threads} {params.extra} {params.index} {input.reads}"
        " | {params.bamsormadup_singularity} bamsormadup {params.tmp_dir} inputformat=sam threads={threads} outputformat=bam level=0 SO=coordinate"
        " | {params.umis_singularity} umis bamtag -"
        " | {params.samtools_singularity} samtools view -b -o {output} - ) &> {log}"
    # wrapper:
    #     "0.38.0/bio/bwa/mem"

rule samtools_index:
    input:
        "bam/{sample}-sort.bam"
    output:
        "bam/{sample}-sort.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}.log"
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools index {input} {output}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/samtools/index"



#fastq1_files = ["fastq/RNA/" + s + "_R1.fastq.gz" for s in config["RNA_Samples"]]
#fastq2_files = ["fastq/RNA/" + s + "_R2.fastq.gz" for s in config["RNA_Samples"]]


rule cutadapt:
    input:
        fastq1 = "fastq/RNA/{sample}_R1_untrimmed.fastq.gz",
        fastq2 = "fastq/RNA/{sample}_R2_untrimmed.fastq.gz"
    output:
        fastq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq/RNA/{sample}_R2.fastq.gz",
        qc = "fastq/RNA/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters_r1 = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapters_r2 = "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        #others = "--minimum-length 2 -q 20"
        others = "--minimum-length 2"
    threads: 8
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cutadaptv2.5-0.simg"
    wrapper:
        "0.38.0/bio/cutadapt/pe"

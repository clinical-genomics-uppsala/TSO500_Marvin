
rule fix_bed_file :
    input:
        bed = "DATA/TST500C_manifest.bed"
    output:
        bed = "bed/ONCOCNV.bed"
    shell:
        "awk 'BEGIN{{ OFS=\"\t\"}}{{ print $1, $2, $3, NR, \"0\", $4 }}' {input.bed} > {output.bed}"


rule Tumor_levels:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        stats = "ref/Control_stats.txt"
    output:
        stats = "stats/{sample}.stats.txt"
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg"
    shell:
        "perl ONCOCNV/ONCOCNV_getCounts.pl getSampleStats "
        "-m Ampli "
        "-c {input.stats} "
        "-s {input.bam} "
        "-o {output.stats}"

rule Target_bed:
    input:
        stats = "ref/Control_stats.txt"
    output:
        bed = "bed/target.bed"
    shell:
        "cat {input.stats} | grep -v start | awk '{{print $1,$2,$3}}' "
        "| sed \"s/ /\t/g\" >{output.bed}"

rule Target_GC:
    input:
        bed = "bed/target.bed",
        reference = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta"
    output:
        GC = "ref/target.GC.txt"
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg"
    shell:
        "perl ONCOCNV/createTargetGC.pl "
        "-bed {input.bed} "
        "-fi {input.reference} "
        "-od stats/ "
        "-of {output.GC}"


rule Calls:
    input:
        tumor_stats = "stats/{sample}.stats.txt",
        control_stats = "ref/Control_stats.Processed.txt"
    output:
        call = "CNV_calls/{sample}.output.txt"
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg"
    shell:
        "cat ONCOCNV/processSamples.R | R --slave "
        "--args {input.tumor_stats} {input.control_stats} {output.call} "
        " cghseg"

rule CNV_event:
    input:
        calls = ["CNV_calls/" + sample_id + ".output.txt" for sample_id in config["DNA_Samples"]]
    output:
        cnv_event = "CNV_calls/cnv_event.txt"
    shell:
        "python src/cnv_event.py "
        "{input.calls} "
        "{output.cnv_event}"

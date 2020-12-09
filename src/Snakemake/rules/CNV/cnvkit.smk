



rule Create_targets:
    input:
        bed = "DATA/TST500C_manifest.bed"
    output:
        bed = "bed/manifest.target.bed"
    threads: 1
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cnvkit.simg"
    shell:
        "cnvkit.py target "
        "--split "
        "{input.bed} "
        "-o {output.bed}"

rule Create_anti_targets:
    input:
        bed = "bed/manifest.target.bed"
    output:
        bed = "bed/manifest.antitarget.bed"
    threads: 1
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cnvkit.simg"
    shell:
        "cnvkit.py antitarget "
        "{input.bed} "
        "-o {output.bed}"



rule Call_cnv:
    input:
        #bams = expand("{tumor_sample}", tumor_sample=config["Tumor_samples"].values()),
        bams = ["DNA_bam/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        cnv_reference = "ref/normal_reference.cnn"
    output:
        #regions = ["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["Tumor_samples"]],
        regions = ["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["DNA_Samples"]],
        #segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["Tumor_samples"]]
        segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]],
    threads: 8
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cnvkit.simg"
    shell:
        "cnvkit.py batch "
        "{input.bams} "
        "-r {input.cnv_reference} "
        "-p {threads} "
        "-d CNV_calls/"

rule Filter_cnv:
    input:
        #segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["Tumor_samples"]],
        segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]],
        purity = "DATA/Pathological_purity_BMS_validation.txt",
        relevant_genes = "DATA/TSO500_relevant_genes.txt",
        ONCOCNV_events = "CNV_calls/cnv_event.txt",
        bed_file = "bed/manifest.target.bed",
        vcf_files = ["Results/DNA/" + sample_id + "/vcf/" + sample_id + "-ensemble.final.vcf.gz" for sample_id in config["DNA_Samples"]]
    output:
        relevant_cnvs = "CNV_results/relevant_cnv.txt",
        cnv_done = "CNV_results/cnv_done.txt"
    #singularity:
    #    "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cnvkit_0.9.7--py_1.sif"
    shell:
        #"source /projects/wp4/nobackup/workspace/jonas_test/CNV_runs/20190909_JA_refernce/CNV/Helena/cnvkit_venv/bin/activate; "
        "python src/report_cnv.py "
        "TSO500 "
        "{input.purity} "
        "{input.relevant_genes} "
        "{input.segments} "
        "{input.ONCOCNV_events} "
        "{input.bed_file} "
        "{output.relevant_cnvs}"


#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./cnvkit.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects --bind /beegfs " --restart-times 2

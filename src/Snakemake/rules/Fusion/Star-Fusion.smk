
#singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg"

rule STAR:
    input:
        fq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fq2 = "fastq/RNA/{sample}_R2.fastq.gz"
    output:
        alignment = "RNA_bam/{sample}Chimeric.out.junction"
    threads: 5
    run:
        import subprocess
        command = "singularity exec -B /projects/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
        #command += "STAR --genomeDir /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/ "
        command += "STAR --genomeDir /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/ "
        command += "--readFilesIn " + input.fq1 + " " + input.fq2 + " "
        command += "--outReadsUnmapped None "
        command += "--twopassMode Basic "
        command += "--readFilesCommand \"gunzip -c\" "
        command += "--outSAMstrandField intronMotif "  # include for potential use with StringTie for assembly
        command += "--outSAMtype BAM SortedByCoordinate "
        command += "--outSAMunmapped Within "
        command += "--chimSegmentMin 12 "  # ** essential to invoke chimeric read detection & reporting **
        command += "--chimJunctionOverhangMin 12 "
        command += "--chimOutJunctionFormat 1 "   # **essential** includes required metadata in Chimeric.junction.out file.
        command += "--alignSJDBoverhangMin 10 "
        command += "--alignMatesGapMax 100000 "   # avoid readthru fusions within 100k
        command += "--alignIntronMax 100000 "
        command += "--alignSJstitchMismatchNmax 5 -1 5 5 "   # settings improved certain chimera detections
        command += "--outSAMattrRGline ID:GRPundef "
        command += "--chimMultimapScoreRange 3 "
        command += "--chimScoreJunctionNonGTAG -4 "
        command += "--chimMultimapNmax 20 "
        command += "--chimNonchimScoreDropMin 10 "
        command += "--peOverlapNbasesMin 12 "
        command += "--peOverlapMMp 0.1 "
        command += "--runThreadN " + str(threads) + " "
        command += "--outFileNamePrefix RNA_bam/" + wildcards.sample
        print(command)
        subprocess.call(command, shell=True)

# rule STAR_Fusion:
#     input:
#         alignment = "STAR2/{sample}Chimeric.out.junction"
#     output:
#         #fusion1 = "Results/RNA/{sample}/Fusions/star-fusion.fusion_predictions.tsv",
#         #fusion2 = "Results/RNA/{sample}/Fusions/star-fusion.fusion_predictions.abridged.tsv"
#         #fusion1 = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.tsv",
#         fusion2 = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.abridged.tsv"
#     threads: 5
#     shell:
#         "singularity exec -B /projects/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
#         "/usr/local/src/STAR-Fusion/STAR-Fusion "
#         #"--genome_lib_dir /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ "
#         "--genome_lib_dir /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ "
#         "-J {input.alignment} "
#         "--output_dir STAR_fusion/{wildcards.sample}/Fusions/ "
#         "--CPU {threads}"


rule STAR_Fusion:
    input:
        alignment = "RNA_bam/{sample}Chimeric.out.junction",
        fq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fq2 = "fastq/RNA/{sample}_R2.fastq.gz"
    output:
        fusion = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
        html = "STAR_fusion/{sample}/Fusions/FusionInspector-inspect/finspector.fusion_inspector_web.html"
    params:
        ref = config["reference"]["STAR_fusion"]
    singularity:
        config["singularity"]["STAR_fusion"]
    threads: 5
    shell:
        #"singularity exec -B /projects/ -B /scratch/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
        "/usr/local/src/STAR-Fusion/STAR-Fusion "
        "--genome_lib_dir {params.ref} "
        "-J {input.alignment} "
        "--output_dir STAR_fusion/{wildcards.sample}/Fusions/ "
        "--CPU {threads} "
        "--left_fq {input.fq1} "
        "--right_fq {input.fq2} "
        "--examine_coding_effect "
        "--FusionInspector inspect"

# rule Star_fusion_validate:
#     input:
#         fusion = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.abridged.tsv",
#         fq1 = "fastq/RNA/{sample}_R1.fastq.gz",
#         fq2 = "fastq/RNA/{sample}_R2.fastq.gz"
#     output:
#         fusion = "FI/{sample}/Fusions/FI/finspector/finspector.FusionInspector.fusions.abridged.tsv"
#     run:
#         import subprocess
#         command = "singularity exec -B /projects/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
#         command += "/usr/local/src/STAR-Fusion/FusionInspector/FusionInspector --fusions STAR_fusion/" + wildcards.sample + "/Fusions/star-fusion.fusion_predictions.abridged.tsv "
#         #command += "--genome_lib /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ "
#         command += "--genome_lib /projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ "
#         command += "--left_fq " + input.fq1 + " "
#         command += "--right_fq " + input.fq2 + " "
#         command += "--output_dir FI/" + wildcards.sample + "/Fusions/FI/finspector "
#         command += "--vis"
#         print(command)
#         subprocess.call(command, shell=True)

rule Copy_to_results:
    input:
        #STAR_fusion1 = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.tsv",
        STAR_fusion2 = "STAR_fusion/{sample}/Fusions/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
        #FI = "FI/{sample}/Fusions/FI/finspector/finspector.FusionInspector.fusions.abridged.tsv"
        html = "STAR_fusion/{sample}/Fusions/FusionInspector-inspect/finspector.fusion_inspector_web.html"
    output:
        #STAR_fusion1 = "Results/RNA/{sample}/Fusions/star-fusion.fusion_predictions.tsv",
        STAR_fusion2 = "Results/RNA/{sample}/Fusions/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
        #FI = "Results/RNA/{sample}/Fusions/finspector.FusionInspector.fusions.abridged.tsv"
        html = "Results/RNA/{sample}/Fusions/Fusion_inspector_web.html"
    shell:
        "cp {input.STAR_fusion2} {output.STAR_fusion2} && "
        "cp {input.html} {output.html}"
        #"cp {input.FI} {output.FI}"



#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time} -w compute01"  -s ./CNV_TSO500.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects  " --cluster-config Config/Slurm/cluster.json --restart-times 2
#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 6 -t 24:00:00"  -s ./Star-Fusion.smk --use-singularity --singularity-args " --bind /projects  "

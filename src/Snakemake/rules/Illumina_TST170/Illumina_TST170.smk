

rule run_TST170:
    output:
        fusions = ["RNA_TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["RNA_TST170/bam_files/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["RNA_TST170/bam_files/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    params:
        runfolder = config["Runfolder"],
        samples = config["RNA_Samples"]
    run:
        import subprocess
        import os
        subprocess.call("mkdir TST170/", shell=True)
        subprocess.call("singularity run -e -B /etc/localtime:/etc/localtime -B " + params.runfolder + ":/data -B /data/illumina/TST170/resources_TST170/genomes:/genomes -B ./TST170:/analysis /projects/wp4/nobackup/workspace/somatic_dev/singularity/docker-oncology.dockerhub.illumina.com_tst170localapp_1.0.0.0-2017-07-28-71e1b6fbab65.sif", shell=True)
        TST170_outfolder = [i for i in os.listdir("TST170/") if "TruSightTumor170_Analysis_" in i][0]
        subprocess.call("mv TST170/" + TST170_outfolder + "/RNA_IntermediateFiles/Alignment/* RNA_TST170/bam_files/", shell=True)
        subprocess.call("rm RNA_TST170/bam_files/*.bed", shell=True)
        subprocess.call("mv TST170/" + TST170_outfolder + "/*.* RNA_TST170/", shell=True)
        for sample in params.samples :
            subprocess.call("mv TST170/" + TST170_outfolder + "/RNA_" + sample + "/* RNA_TST170/RNA_" + sample + "/", shell=True)
        subprocess.call("cp *_samplesheet.csv " + params.runfolder + "SampleSheet.txt", shell=True)

rule TST170_QC_coverage:
    input:
        bam = "RNA_TST170/bam_files/{sample}.bam",
        bed = "DATA/TST500C_manifest.bed"
    output:
        coverage = "Results/RNA/{sample}/{sample}.Housekeeping_gene_coverage.txt"
    run:
        import subprocess
        subprocess.call("python src/RNA_coverage.py " + input.bam + " " + input.bed + " " + output.coverage, shell=True)

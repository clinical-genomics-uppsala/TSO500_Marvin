
rule run_TSO500:
    input:
        sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    output:
         TSO500_done = "TSO500/TSO500_done.txt"
    params:
        runfolder = config["Runfolder"][:-1],
        folder = config["Runfolder"].split("/")[-2]
    shell:
        #"mkdir /scratch/TSO500/ && "
        #"rsync -Prv {params.runfolder} /scratch/ && "
        #"/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ "
        #"--analysisFolder /scratch/TSO500/ --runFolder /scratch/{params.folder}/ --sampleSheet {input.sample_sheet} && "
        "/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ "
        "--analysisFolder TSO500 --runFolder {params.runfolder} --sampleSheet {input.sample_sheet} && "
        #"rsync -Prv /scratch/TSO500/* TSO500/ && "
        #"rm -r /scratch/TSO500/ && "
        #"rm -r /scratch/{params.folder}/ && "
        "touch TSO500/TSO500_done.txt"


localrules: copy_TST170


rule copy_TSO500_RNA:
    input:
        #fusions = ["TSO500/Results/" + s + "/" + s + "_CombinedVariantOutput.tsv" for s in config["RNA_Samples"]],
        TSO500_done = "TSO500_done.txt"
    output:
        fusions = ["Results/RNA/" + s + "/Fusions/" + s + "_CombinedVariantOutput.tsv" for s in config["RNA_Samples"]],
        #bams = ["RNA_TSO500/bam_files/" + s + ".bam" for s in config["RNA_Samples"]],
        #bais = ["RNA_TSO500/bam_files/" + s + ".bam.bai" for s in config["RNA_Samples"]],
        metrics = "Results/RNA/MetricsOutput.tsv"
    params:
        samples = config["RNA_Samples"]
    run:
        import subprocess
        subprocess.call("cp TSO500/Results/MetricsOutput.tsv Results/RNA/", shell=True)
        for sample in params.samples :
            subprocess.call("cp TSO500/Results/" + sample + "/" + sample + "_CombinedVariantOutput.tsv Results/RNA/" + sample + "/Fusions/", shell=True)
            subprocess.call("mv TSO500/Logs_Intermediates/RnaAlignment/" + sample + "/*.bam* RNA_TSO500/bam_files/", shell=True)


# rule copy_IGV:
#     input:
#         done_file = "Arriba_results/{sample}/IGV_done.txt"
#     output:
#         done_file = "Arriba_results/{sample}/IGV_copy_done.txt"
#     shell:
#         "cp Arriba_results/{wildcards.sample}/*.svg Results/RNA/{wildcards.sample}/Fusions/"

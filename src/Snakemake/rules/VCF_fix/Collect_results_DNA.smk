
rule ensemble_filter:
    input:
        vcf = "recall/{sample}.ensemble.vcf.gz"
    output:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.vcf.gz"
    container:
        config["singularity"]["bcftools"]
    shell:
        "python3 src/filter_by_num_callers.py -v {input.vcf} -d | bgzip > {output.vcf} && "
        "tabix {output.vcf}"

rule intron_filter:
    input:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.vcf.gz"
    output:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    params:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf"
    container:
        config["singularity"]["bcftools"]
    shell :
        "python3 src/filter_TSO500_introns.py {input.vcf} &&"
        "bgzip {params.vcf} && "
        "tabix {output.vcf}"

rule AD_filter:
    input:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    output:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vcf"
    container:
        config["singularity"]["bcftools"]
    shell :
        "bcftools filter -O v -o {output.vcf} -e \"FORMAT/AD[0:0]<20\" {input.vcf}"


rule vep:
    input:
        vcf="Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vcf",
        cache=config["configfiles"]["vep"],
        fasta=config["reference"]["ref"],
    output:
        vcf="Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.vcf",
    params:
        # "--everything --check_existing --pick"  #--exclude_null_alleles
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class ",
    log:
        "logs/variantCalling/vep/{sample}.log",
    singularity:
        config["singularity"]["vep"]
    threads: 10
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --refseq --offline --fasta {input.fasta} {params} ) &> {log}"


rule Find_multibp_SNV:
    input:
        vcf="Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.vcf",
        ref=config["reference"]["ref"],
    output:
        vcf=temp("Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.multiplebp.vcf.temp"),
    #singularity:
    #    config["singularity"].get("python_samtools", config["singularity"].get("default", ""))
    script:
        "../../../Multibp_SNV.py"


rule sort_multiplebp_vcf:
    input:
        vcf="Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.multiplebp.vcf.temp",
    output:
        vcf="Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.multiplebp.vcf",
    container:
        config["singularity"]["bcftools"]
    wrapper:
        "0.72.0/bio/bcftools/sort"


rule ffpe_filter:
    input:
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vcf",
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai"
    params:
        vcf_ffpe_temp = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.ffpe.temp.vcf",
        vcf_ffpe = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.ffpe.tsv"
    output:
        #gvcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vcf.gz",
        gvcf_ffpe = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.ffpe.tsv.gz"
    container:
        config["singularity"]["bcftools"]
    shell:
        #"module load oracle-jdk-1.8/1.8.0_162 && "
        "java -jar SOBDetector/SOBDetector_v1.0.1.jar --input-type VCF --input-variants {input.vcf} --input-bam {input.bam} --output-variants {params.vcf_ffpe_temp} && "
        "python src/Add_FFPE_column_to_vcf.py {params.vcf_ffpe_temp} {params.vcf_ffpe} && "
        "rm {params.vcf_ffpe_temp} && "
        "bgzip {params.vcf_ffpe} && "
        "tabix {output.gvcf_ffpe}"
        #"bgzip {input.vcf} && "
        #"tabix {output.gvcf}"


rule copy_mv_TS0500:
    input:
        TSO500_done = "TSO500/TSO500_done.txt"
    output:
        biomarker = ["Results/DNA/" + s + "/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
        metrics = "Results/DNA/MetricsReport.tsv"
    params:
        samples = config["DNA_Samples"]
    run:
        import subprocess
        subprocess.call("mkdir DNA_TSO500/", shell=True)
        subprocess.call("mkdir DNA_TSO500/Fastq/", shell=True)
        for sample in params.samples :
            subprocess.call("cp TSO500/Results/" + sample + "_BiomarkerReport.txt Results/DNA/" + sample + "/", shell=True)
            subprocess.call("mkdir DNA_TSO500/Fastq/" + sample + "/", shell=True)
            subprocess.call("mv TSO500/Logs_Intermediates/FastqGeneration/" + sample + "/* DNA_TSO500/Fastq/" + sample + "/", shell=True)
        subprocess.call("cp TSO500/Results/MetricsReport.tsv Results/DNA/", shell=True)
        subprocess.call("mv TSO500/Logs_Intermediates/Tmb/ DNA_TSO500/Tmb/", shell=True)
        subprocess.call("mv TSO500/Logs_Intermediates/Msi/ DNA_TSO500/Msi/", shell=True)
        subprocess.call("mv TSO500/Logs_Intermediates/StitchedReads/ DNA_TSO500/StitchedReads/", shell=True)
        subprocess.call("mv TSO500/Logs_Intermediates/RunQc/ DNA_TSO500/RunQc/", shell=True)
        subprocess.call("mv TSO500/Logs_Intermediates/VariantCaller/ DNA_TSO500/VariantCaller/", shell=True)


rule copy_CNV:
    input:
        cnv = "CNV_results/relevant_cnv.txt",
        cnv_done = "CNV_results/cnv_done.txt"
    output:
        cnv_png = ["Results/DNA/" + s + "/CNV/" + s + "-ready.png" for s in config["DNA_Samples"]]
    params:
        DNA_samples = expand(config["DNA_Samples"])
    run:
        import subprocess
        for sample in params.DNA_samples :
            #subprocess.call("cp CNV_results/relevant_cnv.txt Results/DNA/" + sample + "/CNV/", shell=True)
            subprocess.call("grep \"sample_path\" CNV_results/relevant_cnv.txt > Results/DNA/" + sample + "/CNV/relevant_cnv.txt", shell=True)
            subprocess.call("grep \"" + sample + "-ready.cns\" CNV_results/relevant_cnv.txt >> Results/DNA/" + sample + "/CNV/relevant_cnv.txt", shell=True)
            subprocess.call("cp CNV_results/" + sample + "-ready*.png " + "Results/DNA/" + sample + "/CNV/", shell=True)

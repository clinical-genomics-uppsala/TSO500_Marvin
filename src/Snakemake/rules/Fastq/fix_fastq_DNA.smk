S_dna = []
for s in config["DNA_Samples"].values() :
    S_dna.append(s)
fastq1_files = ["fastq_temp/DNA/" + s + "_" + i + "_R1_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]
fastq2_files =  ["fastq_temp/DNA/" + s + "_" + i + "_R2_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]


rule fix_fastq_bash_DNA:
    input:
        fastq1 = fastq1_files,
        fastq2 = fastq2_files
    output:
        bash_scripts_dna_R1 = ["fastq_temp/DNA/" + s + "_R1.fix_fastq.sh" for s in config["DNA_Samples"]],
        bash_scripts_dna_R2 = ["fastq_temp/DNA/" + s + "_R2.fix_fastq.sh" for s in config["DNA_Samples"]]
    params:
        DNA_samples = [s for s in config["DNA_Samples"]]
    run:
        import subprocess
        subprocess.call("mkdir fastq",shell=True)
        i = 0
        for sample in params.DNA_samples :
            bs = open("fastq_temp/DNA/" + sample + "_R1.fix_fastq.sh", "w")
            bs.write("for s in " + S_dna[i] + "," + sample + "; do\n")
            bs.write("\tIFS=\",\";\n")
            bs.write("\tset -- $s;\n")
            bs.write("\tsample_number=$1;\n")
            bs.write("\tsample=$2\n")
            bs.write("\t\tfor r in R1; do\n")
            bs.write("\t\t\techo \"zcat fastq_temp/DNA/\"$sample\"_\"$sample_number\"_\"$r\"* | awk '{if(/^@/){split(\$0,a,\\\":\\\");print(a[1]\\\":\\\"a[2]\\\":\\\"a[3]\\\":\\\"a[4]\\\":\\\"a[5]\\\":\\\"a[6]\\\":\\\"a[7]\\\":UMI_\\\"gsub(\\\"+\\\",\\\"\\\",a[8])\\\":\\\"a[9]\\\":\\\"a[10]\\\":\\\"a[11])}else{print(\$0)}}' | gzip > fastq/DNA/\"$sample\"_\"$r\".fastq.gz \";\n")
            bs.write("\t\tdone  | bash -\n")
            bs.write("done\n")
            #bs.write("sleep 7100\n")
            bs.close()
            subprocess.call("chmod 774 fastq_temp/DNA/" + sample + "_R1.fix_fastq.sh", shell=True)
            i += 1
        i = 0
        for sample in params.DNA_samples :
            bs = open("fastq_temp/DNA/" + sample + "_R2.fix_fastq.sh", "w")
            bs.write("for s in " + S_dna[i] + "," + sample + "; do\n")
            bs.write("\tIFS=\",\";\n")
            bs.write("\tset -- $s;\n")
            bs.write("\tsample_number=$1;\n")
            bs.write("\tsample=$2\n")
            bs.write("\t\tfor r in R2; do\n")
            bs.write("\t\t\techo \"zcat fastq_temp/DNA/\"$sample\"_\"$sample_number\"_\"$r\"* | awk '{if(/^@/){split(\$0,a,\\\":\\\");print(a[1]\\\":\\\"a[2]\\\":\\\"a[3]\\\":\\\"a[4]\\\":\\\"a[5]\\\":\\\"a[6]\\\":\\\"a[7]\\\":UMI_\\\"gsub(\\\"+\\\",\\\"\\\",a[8])\\\":\\\"a[9]\\\":\\\"a[10]\\\":\\\"a[11])}else{print(\$0)}}' | gzip > fastq/DNA/\"$sample\"_\"$r\".fastq.gz \";\n")
            bs.write("\t\tdone  | bash -\n")
            bs.write("done\n")
            #bs.write("sleep 7100\n")
            bs.close()
            subprocess.call("chmod 774 fastq_temp/DNA/" + sample + "_R2.fix_fastq.sh", shell=True)
            i += 1



rule fix_fastq_run_DNA_R1:
    input:
        bash_scripts_DNA_R1 = "fastq_temp/DNA/{sample}_R1.fix_fastq.sh",
        fastq1 = fastq1_files,
        fastq2 = fastq2_files,
    output:
        merged_fastq_R1_DNA = "fastq/DNA/{sample}_R1.fastq.gz"
    shell:
        "{input.bash_scripts_DNA_R1}"

rule fix_fastq_run_DNA_R2:
    input:
        bash_scripts_DNA_R2 = "fastq_temp/DNA/{sample}_R2.fix_fastq.sh",
        fastq1 = fastq1_files,
        fastq2 = fastq2_files,
    output:
        merged_fastq_R2_DNA = "fastq/DNA/{sample}_R2.fastq.gz"
    shell:
        "{input.bash_scripts_DNA_R2}"

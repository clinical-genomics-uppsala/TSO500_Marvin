
localrules: Arriba_HC, Arriba_IGV_bat

rule STAR_arrbia:
    input:
        fastq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq/RNA/{sample}_R2.fastq.gz"
    output:
        bams = "STAR/{sample}Aligned.sortedByCoord.out.bam",
        bais = "STAR/{sample}Aligned.sortedByCoord.out.bam.bai",
        junctions = "STAR/{sample}SJ.out.tab"
    threads: 5
    run:
        import subprocess
        #command = "singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg "
        command = "singularity exec -B /projects/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg "
        command += "STAR "
    	command += "--runThreadN " + str(threads) + " "
    	command += "--genomeDir /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/STAR_index_hs37d5_GENCODE19/ "
        command += "--genomeLoad NoSharedMemory "
    	command += "--readFilesIn " + input.fastq1 + " " + input.fastq2 + " "
        command += "--readFilesCommand zcat "
        command += "--outSAMtype BAM SortedByCoordinate "
        command += "--outSAMunmapped Within "
    	command += "--outFilterMultimapNmax 1 "
        command += "--outFilterMismatchNmax 3 "
    	command += "--chimSegmentMin 10 "
        command += "--chimOutType WithinBAM SoftClip "
        command += "--chimJunctionOverhangMin 10 "
        command += "--chimScoreMin 1 "
        command += "--chimScoreDropMax 30 "
        command += "--chimScoreJunctionNonGTAG 0 "
        command += "--chimScoreSeparation 1 "
        command += "--alignSJstitchMismatchNmax 5 -1 5 5 "
        command += "--chimSegmentReadGapMax 3 "
        command += "--outFileNamePrefix STAR/" + wildcards.sample
        print(command)
        subprocess.call(command, shell=True)
        subprocess.call("samtools index " + output.bams, shell=True)


rule Arriba:
    input:
        bams = "STAR/{sample}Aligned.sortedByCoord.out.bam"
    output:
        fusions1 = "Arriba_results/{sample}.fusions.tsv",
        fusions2 = "Arriba_results/{sample}.fusions.discarded.tsv"
    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg"
    shell:
        "/arriba_v1.1.0/arriba "
    	"-x {input.bams} "
    	"-o {output.fusions1} "
        "-O {output.fusions2} "
    	"-a /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/hs37d5.fa "
        "-g /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/GENCODE19.gtf "
        "-b /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/blacklist_hg19_hs37d5_GRCh37_2018-11-04.tsv "
    	"-T "
        "-P "

rule Arriba_HC:
    input:
        fusions = "Arriba_results/{sample}.fusions.tsv",
        refseq = "DATA/refseq_full_hg19.txt"
    output:
        #fusions = "Results/RNA/{sample}/Fusions/{sample}.Arriba.HighConfidence.fusions.tsv"
        fusions = "Results/RNA/{sample}/Fusions/{sample}.Arriba.fusions.tsv"
    shell:
        #"head -n 1 {input.fusions} > {output.fusions} && "
        #"grep 'high' {input.fusions} >> {output.fusions} || true && "
        #"python src/Add_fusion_exon_name.py {input.refseq} {output.fusions}"
        "python src/Add_fusion_exon_name.py {input.refseq} {input.fusions} {output.fusions}"


rule Arriba_image:
    input:
        #fusion = "Results/RNA/{sample}/Fusions/{sample}.Arriba.HighConfidence.fusions.tsv",
        fusion = "Results/RNA/{sample}/Fusions/{sample}.Arriba.fusions.tsv",
        bam = "STAR/{sample}Aligned.sortedByCoord.out.bam",
        bai = "STAR/{sample}Aligned.sortedByCoord.out.bam.bai"
    output:
        image = "Results/RNA/{sample}/Fusions/{sample}.Arriba.fusions.pdf"
    params:
        image_out_path = "Results/RNA/{sample}/Fusions/"
    run:
        import subprocess
        command = "singularity exec "
        command += "-B " + params.image_out_path + ":/output "
        command += "-B /projects/wp4/nobackup/workspace/jonas_test/Arriba/references:/references:ro "
        command += "-B " + input.fusion + ":/fusions.tsv:ro "
        command += "-B " + input.bam + ":/Aligned.sortedByCoord.out.bam:ro "
        command += "-B " + input.bai + ":/Aligned.sortedByCoord.out.bam.bai:ro "
        command += "/projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg "
        command += "draw_fusions.sh"
        print(command)
        subprocess.call(command, shell=True)
        subprocess.call("mv " + params.image_out_path + "fusions.pdf " + output.image, shell=True)

#rule Arriba_IGV_bat:
#    input:
#        fusion = "Results/RNA/{sample}/{sample}.Arriba.HighConfidence.fusions.tsv",
#        bam = "STAR/{sample}Aligned.sortedByCoord.out.bam"
#    output:
#        bat = "Arriba_results/{sample}_IGV.bat"
#    run:
#        sample = output.bat.split("Arriba_results/")[1].split("_IGV")[0]
#        outfile = open(output.bat, "w")
#        outfile.write("new\n")
#        outfile.write("genome DATA/igv/genomes/hg19.genome\n")
#        outfile.write("load " + input.bam + "\n")
#        outfile.write("snapshotDirectory Arriba_results/" + sample + "/\n")
#        header = True
#        infile = open(input.fusion)
#        for line in infile :
#            if header:
#                header = False
#                continue
#            lline = line.strip().split("\t")
#            pos1 = "chr" + lline[4]
#            pos2 = "chr" + lline[5]
#            gene1 = lline[0]
#            gene2 = lline[1]
#            outfile.write("goto " + pos1 + " " + pos2 + "\n")
#            outfile.write("sort base\n")
#            outfile.write("squish\n")
#            outfile.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
#            outfile.write("snapshot " + gene1 + "_" + gene2 + "_" + pos1.split(":")[0] + "_" + pos1.split(":")[1] + "_" + pos2.split(":")[0] + "_" + pos2.split(":")[1] + ".svg\n")
#        outfile.write("exit\n")
#        infile.close()
#        outfile.close()

#rule Arriba_IGV_run:
#    input:
#        bat = "Arriba_results/{sample}_IGV.bat"
#    output:
#        done_file = "Arriba_results/{sample}/IGV_done.txt"
#    singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/igv-2.4.10-0.simg"
#    shell:
#        "xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} && touch {output.done_file}"

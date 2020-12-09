
rule varscan:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        ref = config["reference"]["ref"],
        bed = config["bed"]["bedfile"]
    output:
        temp("varscan/{sample}.varscan.vcf")
    params:
        samtools_singularity = config["singularity"]["execute"] + config["singularity"]["samtools"],
        varscan_singularity = config["singularity"]["execute"] + config["singularity"]["varscan"],
        mpileup = "-d 1000 -L 1000",
        varscan = "--min-coverage 5 --p-value 0.98 --strand-filter 1 --min-var-freq 0.01 --output-vcf --variants"
    log:
        "logs/variantCalling/varscan/{sample}.log"
    shell:
        "({params.samtools_singularity} samtools mpileup -f {input.ref} {params.mpileup} -l {input.bed} {input.bam} |"
        " grep -v -P '\t0\t\t$' |"
        " {params.varscan_singularity} java -jar /usr/local/share/varscan-2.4.3-0/VarScan.jar mpileup2cns {params.varscan} |" #" --vcf-sample-list sample_list.txt ""
        #" | /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/py -x 'bcbio.variation.vcfutils.add_contig_to_header(x, "/data/ref_genomes/bcbio-nextgen/sam/hg19.with.mt.fasta")'
        #" | /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/py -x 'bcbio.variation.varscan.fix_varscan_output(x)' | "
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' | "
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $5) }} {{print}}' | "
        " sed 's/:\([0-9]\+\)\.\([0-9]\+\)%:/:0.\\1\\2:/' | " #Exchange alleles frequences from ex 10.5% to 0.105
        " sed 's/:\([0-9]\+\)%:/:0.\\1:/' | " #Exchange alleles frequences from ex 10% to 0.10
        " sed 's/:100%:/:1.00:/' " #Exchange alleles frequences from 100% to 1.00

        #" sed 's/FREQ/AF/' " #Exchange FORMAT/FREQ to FORMAT/AF
        " > {output}) 2> {log}"
        #" ifne vcfuniqalleles > /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmpFclrjK/20-1618-chr7_0_41726120.vcf"

rule sortVarscan:
    input:
        "varscan/{sample}.varscan.vcf"
    output:
        temp("varscan/{sample}.varscan.fixAF.vcf")
    singularity:
        config["singularity"]["bcftools"]
    log:
        "logs/variantCalling/varscan/{sample}.sort.log"
    shell:
        "(bgzip {input} && tabix {input}.gz && "
        "bcftools sort -o {output} -O v {input}.gz) &> {log}"


#gunzip -c /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/gemini/20-1578-varscan.vcf.gz | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -
#s - | vt normalize -n -r /data/ref_genomes/bcbio-nextgen/sam/hg19.with.mt.fasta - | awk '{ gsub("./-65", "./."); print $0 }'| bcftools view -f 'PASS,.' | sed -e 's/Number=A/Number=1/g' | bgzip -c
# > /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmp1T8c8T/20-1578-varscan-decompose.vcf.gz
#[2020-08-24T08:51Z] compute04: /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/tabix -f -p vcf /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmplmZICZ/20-1578-varscan-
#decompose.vcf.gz

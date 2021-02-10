
include: "../rules/Fastq/demultiplex.smk"

if config["DNA_Samples"] != "No DNA" :
    include: "../rules/Illumina_TSO500/TSO500_SampleSheet.smk"
    include: "../rules/Illumina_TSO500/Illumina_TSO500.smk"
    include: "../rules/Fastq/fix_fastq_DNA.smk"
    #include: "../rules/Bcbio/Bcbio.smk"
    include: "../rules/CNV/ONCOCNV.smk"
    include: "../rules/CNV/cnvkit.smk"
    include: "../rules/QC/check_coverage.smk"
    include: "../rules/VCF_fix/Collect_results_DNA.smk" #Change folder!
    #include: "../rules/Mutect2/Mutect2.smk"

    include: "../rules/Alignment/bwa-mem.smk"
    include: "../rules/Alignment/fgbio.smk"
    include: "../rules/Alignment/GATK.smk"
    include: "../rules/SNV/freebayes.smk"
    include: "../rules/SNV/mutect2.smk"
    include: "../rules/SNV/vardict_T.smk"
    include: "../rules/SNV/varscan.smk"
    include: "../rules/VCF_fix/fix_AF_all_callers.smk"
    include: "../rules/VCF_fix/normalize.smk"
    include: "../rules/VCF_fix/recall.smk"
    include: "../rules/QC/samtools-picard-stats.smk"
    include: "../rules/QC/multiqc.smk"
    include: "../rules/QC/cartool.smk"
    include: "../rules/QC/fastqc.smk"


if config["RNA_Samples"] != "No RNA" :
    include: "../rules/Illumina_TST170/TST170_SampleSheet.smk"
    include: "../rules/Illumina_TST170/Illumina_TST170.smk"
    include: "../rules/Fastq/fix_fastq_RNA.smk"
    include: "../rules/Adapter_trimming/Adapter_trimming.smk"
    include: "../rules/Fusion/Arriba.smk"
    include: "../rules/Fusion/Imbalance.smk"
    #include: "../rules/Fusion/exon_splicing.smk"
    include: "../rules/Fusion/Exon_skipping.smk"
    include: "../rules/Fusion/Star-Fusion.smk"
    include: "../rules/VCF_fix/Collect_results_RNA.smk" #Change folder!
    include: "../rules/QC/RSeQC.smk"

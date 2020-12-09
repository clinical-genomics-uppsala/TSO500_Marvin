
#include: "../rules/Fastq/demultiplex.smk"
include: "../rules/Illumina_TSO500/TSO500_SampleSheet_v2.smk"
include: "../rules/Illumina_TSO500/Illumina_TSO500_v2.smk"

if config["DNA_Samples"] != "No DNA" :
    #include: "../rules/Illumina_TSO500/TSO500_SampleSheet_v2.smk"
    #include: "../rules/Illumina_TSO500/Illumina_TSO500_v2.smk"
    #include: "../rules/Fastq/fix_fastq_DNA.smk"
    #include: "../rules/Bcbio/Bcbio.smk"
    #include: "../rules/CNV/ONCOCNV.smk"
    #include: "../rules/CNV/cnvkit.smk"
    #include: "../rules/QC/check_coverage.smk"
    include: "../rules/VCF_fix/Collect_results_DNA_v2.smk"
    #include: "../rules/Mutect2/Mutect2.smk"

if config["RNA_Samples"] != "No RNA" :
    #include: "../rules/Illumina_TST170/TST170_SampleSheet.smk"
    #include: "../rules/Fastq/fix_fastq_RNA.smk"
    #include: "../rules/Adapter_trimming/Adapter_trimming.smk"
    #include: "../rules/Fusion/Arriba.smk"
    #include: "../rules/Fusion/Imbalance_v2.smk"
    #include: "../rules/Fusion/exon_splicing.smk"
    #include: "../rules/Illumina_TST170/Illumina_TST170.smk"
    include: "../rules/Collect_results/Collect_results_RNA_v2.smk"
    #include: "../rules/QC/RSeQC_v2.smk"
    #include: "../rules/Fusion/Star-Fusion.smk"

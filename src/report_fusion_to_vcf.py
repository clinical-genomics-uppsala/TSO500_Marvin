
#import glob
#import os
import sys
#import subprocess
from datetime import date
from subprocess import check_output


#fusion_file_name = snakemake.input.variants
#fusion_vcf = open(snakemake.ouput.vcf, "w")

fusion_file_name = sys.argv[1]
in_fastq_ref = sys.argv[2]
fusion_vcf = open(sys.argv[3], "w")


#fusion_file_name = "Results/RNA/R21-475/Fusions/R21-475_HighConfidenceVariants.csv"
#fusion_vcf = open("Results/RNA/R21-475/Fusions/R21-475_HighConfidenceVariants.vcf", "w")
fusion_file = open(fusion_file_name)


def write_vcf_header(sample_name, vcf_file):
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##fileDate=%s\n" % str(date.today()))
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_file.write("##INFO=<ID=GENE_NAME1,Number=1,Type=String,Description=\"Fusion gene\">\n")
    vcf_file.write("##INFO=<ID=GENE_NAME2,Number=1,Type=String,Description=\"Fusion partner\">\n")
    vcf_file.write("##INFO=<ID=Fusion_support,Number=1,Type=Integer,Description=\"Number of fusion supporting reads (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Gene1_ref_reads,Number=1,Type=String,Description=\"Number of reads in gene1 not supporting the fusion (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Gene2_ref_reads,Number=1,Type=String,Description=\"Number of reads in gene2 not supporting the fusion (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Breakpoint1,Number=1,Type=String,Description=\"Breakpoint in gene1\">\n")
    vcf_file.write("##INFO=<ID=Breakpoint2,Number=1,Type=String,Description=\"Breakpoint in gene2\">\n")
    vcf_file.write("##INFO=<ID=GENE_NAME,Number=1,Type=String,Description=\"Splice gene\">\n")
    vcf_file.write("##INFO=<ID=Splice_support,Number=1,Type=Integer,Description=\"Number of splice supporting reads (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Splice1_ref_reads,Number=1,Type=String,Description=\"Number of reads in splice cite 1 not supporting the splicing (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Splice2_ref_reads,Number=1,Type=String,Description=\"Number of reads in splice cite 2 not supporting the splicing (including duplicates)\">\n")
    vcf_file.write("##ALT=<ID=Fusion,Description=\"Fusion\">\n")
    vcf_file.write("##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Number of splice supporting reads\">\n")
    vcf_file.write("##FORMAT=<ID=SREF,Number=1,Type=Integer,Description=\"Number of reads supporting the reference\">\n")
    vcf_file.write("##FORMAT=<ID=FS,Number=1,Type=Integer,Description=\"Number of fusion supporting reads\">\n")
    vcf_file.write("##FORMAT=<ID=F1REF,Number=1,Type=Integer,Description=\"Number of reads in gene1 supporting the reference\">\n")
    vcf_file.write("##FORMAT=<ID=F2REF,Number=1,Type=Integer,Description=\"Number of reads in gene2 supporting the reference\">\n")
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)


sample_name = fusion_file_name.split("/")[-1].split("_HighConfidenceVariants")[0]
write_vcf_header(sample_name, fusion_vcf)

'''Read fusion file'''
Fusions = []
state = "Header"
for line in fusion_file :
    if line.find("Gene Fusion") != -1 :
        state = "Fusions"
    elif state == "Fusions" and (line.find("No results found") != -1 or line.strip() == "") :
        state = "Inbetween"
    elif line.find("Gene,Affected Exon(s)") != -1 :
        state = "Splice"
    elif state == "Splice" and (line.find("No results found") != -1 or line.strip() == "") :
        state = "Done"
    else :
        columns = line.strip().split(",")
        if state == "Fusions" :
            gene1 = columns[0].split("/")[0]
            gene2 = columns[0].split("/")[1]
            bp1 = columns[1]
            bp2 = columns[2]
            fusion_support = columns[3]
            gene1_ref_reads = columns[4]
            gene2_ref_reads = columns[5]
            score = columns[6]
            chrom = bp1.split(":")[0]
            start_pos = bp1.split(":")[1]

            fasta_pos = int(start_pos) + 1
            command = "samtools faidx " + in_fastq_ref + " " + chrom + ":" + str(fasta_pos) + "-" + str(fasta_pos)
            ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
            ref = ref_bp.upper()
            alt = "<BND>"

            id = "."
            qual = score
            filter = "PASS"
            info = "SVTYPE=BND;GENE_NAME1=%s;GENE_NAME2=%s;Fusion_support=%s;Gene1_ref_reads=%s;Gene2_ref_reads=%s;Breakpoint1=%s;Breakpoint2=%s" % (
                gene1, gene2, fusion_support, gene1_ref_reads, gene2_ref_reads, bp1, bp2
            )
            format = "FS:G1REF:G2REF"
            data = "%s:%s:%s" % (fusion_support, gene1_ref_reads, gene2_ref_reads)
            out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos, id, ref, alt, qual, filter, info, format, data)
            fusion_vcf.write(out_line)
        elif state == "Splice" :
            print(columns)
            gene = columns[0].split("/")[0]
            exons = columns[1]
            transcript = columns[2]
            bp1 = columns[3]
            bp2 = columns[4]
            splice_support = columns[5]
            reference_reads = columns[6]
            score = columns[7]
            chrom = bp1.split(":")[0]
            start_pos = bp1.split(":")[1]
            fasta_pos = int(start_pos) + 1
            command = "samtools faidx " + in_fastq_ref + " " + chrom + ":" + str(fasta_pos) + "-" + str(fasta_pos)
            ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
            ref = ref_bp.upper()
            alt = "<DEL>"
            id = "."
            qual = score
            filter = "PASS"
            info = "SVTYPE=RNAExonVariant;GENE_NAME=%s;Splice_support=%s;reference_reads=%s;Breakpoint1=%s;Breakpoint2=%s" % (
                gene, splice_support, reference_reads, bp1, bp2
            )
            format = "SS:SREF"
            data = "%s:%s" % (splice_support, reference_reads)
            out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos, id, ref, alt, qual, filter, info, format, data)
            fusion_vcf.write(out_line)

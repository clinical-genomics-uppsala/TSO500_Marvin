
#import glob
#import os
import sys
#import subprocess
from datetime import date
from subprocess import check_output


#fusion_file_name = snakemake.input.variants
#fusion_vcf = open(snakemake.ouput.vcf, "w")

fusion_file_name1 = sys.argv[1]
fusion_file2 = open(sys.argv[2])
in_fastq_ref = sys.argv[3]
fusion_vcf = open(sys.argv[4], "w")


#fusion_file_name = "Results/RNA/R21-475/Fusions/R21-475_HighConfidenceVariants.csv"
#fusion_vcf = open("Results/RNA/R21-475/Fusions/R21-475_HighConfidenceVariants.vcf", "w")
fusion_file1 = open(fusion_file_name1)


def write_vcf_header(sample_name, vcf_file):
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##fileDate=%s\n" % str(date.today()))
    vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_file.write("##INFO=<ID=Fusion_support,Number=1,Type=Integer,Description=\"Number of fusion supporting reads (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Gene1_ref_reads,Number=1,Type=String,Description=\"Number of reads in gene1 not supporting the fusion (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Gene2_ref_reads,Number=1,Type=String,Description=\"Number of reads in gene2 not supporting the fusion (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Breakpoint1,Number=1,Type=String,Description=\"Breakpoint in gene1\">\n")
    vcf_file.write("##INFO=<ID=Breakpoint2,Number=1,Type=String,Description=\"Breakpoint in gene2\">\n")
    vcf_file.write("##INFO=<ID=GENE_NAME,Number=1,Type=String,Description=\"Splice or fusion gene\">\n")
    vcf_file.write("##INFO=<ID=Splice_support,Number=1,Type=Integer,Description=\"Number of splice supporting reads (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Splice1_ref_reads,Number=1,Type=String,Description=\"Number of reads in splice cite 1 not supporting the splicing (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=Splice2_ref_reads,Number=1,Type=String,Description=\"Number of reads in splice cite 2 not supporting the splicing (including duplicates)\">\n")
    vcf_file.write("##INFO=<ID=END,Number=1,Type=String,Description=\"End of splicing\">\n")
    vcf_file.write("##ALT=<ID=Fusion,Description=\"Fusion breakpoint and direction\">\n")
    vcf_file.write("##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Number of splice supporting reads\">\n")
    vcf_file.write("##FORMAT=<ID=SREF,Number=1,Type=Integer,Description=\"Number of reads supporting the reference\">\n")
    vcf_file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_file.write("##FORMAT=<ID=FS,Number=1,Type=Integer,Description=\"Number of fusion supporting reads\">\n")
    vcf_file.write("##FORMAT=<ID=F1REF,Number=1,Type=Integer,Description=\"Number of reads in gene1 supporting the reference\">\n")
    vcf_file.write("##FORMAT=<ID=F2REF,Number=1,Type=Integer,Description=\"Number of reads in gene2 supporting the reference\">\n")
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)


sample_name = fusion_file_name1.split("/")[-1].split("_HighConfidenceVariants")[0]
write_vcf_header(sample_name, fusion_vcf)

'''Read fusion file1'''
Fusions = []
state = "Header"
for line in fusion_file1 :
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
            pass
        elif state == "Splice" :
            print(columns)
            gene = columns[0].split("/")[0]
            exons = columns[1]
            transcript = columns[2]
            bp1_list = columns[3].split(":")
            bp2_list = columns[4].split(":")
            bp1 = bp1_list[0] + ":" + str(int(bp1_list[1]))
            bp2 = bp2_list[0] + ":" + str(int(bp2_list[1]))
            splice_support = columns[5]
            reference_reads = columns[6]
            score = columns[7]
            chrom = bp1.split(":")[0]
            start_pos1 = bp1.split(":")[1]
            start_pos2 = bp2.split(":")[1]
            fasta_pos1 = int(start_pos1) - 2
            fasta_pos2 = int(start_pos2) + 2
            command = "samtools faidx " + in_fastq_ref + " " + chrom + ":" + str(fasta_pos1) + "-" + str(fasta_pos1)
            ref_bp1 = check_output(command, shell=True).decode("utf-8").split("\n")[1]
            command = "samtools faidx " + in_fastq_ref + " " + chrom + ":" + str(fasta_pos2) + "-" + str(fasta_pos2)
            ref_bp2 = check_output(command, shell=True).decode("utf-8").split("\n")[1]
            ref1 = ref_bp1.upper()
            ref2 = ref_bp2.upper()
            alt1 = ref1 + "[" + bp2 + "["
            alt2 = "]" + bp1 + "]" + ref2
            id = gene
            qual = score
            filter = "PASS"
            info = "SVTYPE=Fusion;GENE_NAME=%s;Splice_support=%s;reference_reads=%s;END=%s" % (
                gene, splice_support, reference_reads, bp2.split(":")[1]
            )
            format = "SS:SREF:GT:GQ"
            data = "%s:%s:./.:." % (splice_support, reference_reads)
            out_line1 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos1, id, ref1, alt1, qual, filter, info, format, data)
            fusion_vcf.write(out_line1)
            out_line2 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos2, id, ref2, alt2, qual, filter, info, format, data)
            fusion_vcf.write(out_line2)

'''Read fusion file2'''
Fusions = []
header = True
for line in fusion_file2 :
    if header:
        if line.find("Caller,Gene A,Gene B") != -1 :
            header = False
        continue
    columns = line.strip().split(",")
    keepFusion = columns[19]
    if keepFusion != "True" :
        continue
    print(columns)
    gene1 = columns[1]
    gene2 = columns[2]

    bp1_list = columns[3].split(":")
    bp2_list = columns[4].split(":")
    bp1 = bp1_list[0] + ":" + str(int(bp1_list[1]) + 1)
    bp2 = bp2_list[0] + ":" + str(int(bp2_list[1]) + 1)
    chrom1 = bp1.split(":")[0]
    chrom2 = bp2.split(":")[0]
    start_pos1 = bp1.split(":")[1]
    start_pos2 = bp2.split(":")[1]
    command = "samtools faidx " + in_fastq_ref + " " + chrom1 + ":" + str(start_pos1) + "-" + str(start_pos1)
    ref_bp1 = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    command = "samtools faidx " + in_fastq_ref + " " + chrom2 + ":" + str(start_pos2) + "-" + str(start_pos2)
    ref_bp2 = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    ref1 = ref_bp1.upper()
    ref2 = ref_bp2.upper()

    score = columns[5]
    fusion_support = str(int(columns[13]) + int(columns[14]))
    gene1_ref_reads = str(int(columns[9]) + int(columns[10]))
    gene2_ref_reads = str(int(columns[11]) + int(columns[12]))
    contig = columns[16]
    len_gene1_contig = int(columns[17])
    len_gene2_contig = int(columns[18])
    qual = score
    filter = "PASS"

    fasta_pos = int(start_pos1)+1
    command = "samtools faidx " + in_fastq_ref + " " + chrom1 + ":" + str(fasta_pos-19) + "-" + str(fasta_pos)
    ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    ref1_plus1 = ref_bp.upper()[1:-1]
    fasta_pos = int(start_pos1)+1
    command = "samtools faidx " + in_fastq_ref + " " + chrom1 + ":" + str(fasta_pos) + "-" + str(fasta_pos+19)
    ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    ref1_plus2 = ref_bp.upper()[1:-1]
    ref1_minus1 = ""
    ref1_minus2 = ""
    dna_transalation = {"A" : "T", "C" : "G", "T" : "A", "G" : "C", "N" : "N"}
    for bp in ref1_plus1:
        ref1_minus1 = dna_transalation[bp] + ref1_minus1
    for bp in ref1_plus2:
        ref1_minus2 = dna_transalation[bp] + ref1_minus2

    fasta_pos = int(start_pos2)+1
    command = "samtools faidx " + in_fastq_ref + " " + chrom2 + ":" + str(fasta_pos-19) + "-" + str(fasta_pos)
    ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    ref2_plus1 = ref_bp.upper()[1:-1]
    fasta_pos = int(start_pos2)+1
    command = "samtools faidx " + in_fastq_ref + " " + chrom2 + ":" + str(fasta_pos) + "-" + str(fasta_pos+19)
    ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
    ref2_plus2 = ref_bp.upper()[1:-1]
    ref2_minus1 = ""
    ref2_minus2 = ""
    for bp in ref2_plus1:
        ref2_minus1 = dna_transalation[bp] + ref2_minus1
    for bp in ref2_plus2:
        ref2_minus2 = dna_transalation[bp] + ref2_minus2

    g1_pos = ""
    g1_direction = ""
    g2_pos = ""
    g2_direction = ""
    if contig.find(ref1_plus1) != -1 or contig.find(ref1_plus2) != -1 :
        g1_direction = "plus"
        g1_pos = max(contig.find(ref1_plus1), contig.find(ref1_plus2))
    elif contig.find(ref1_minus1) != -1 or contig.find(ref1_minus2) != -1 :
        g1_direction = "minus"
        g1_pos = max(contig.find(ref1_minus1), contig.find(ref1_minus2))
    else:
        g1_direction = "error"

    if contig.find(ref2_plus1) != -1 or contig.find(ref2_plus2) != -1 :
        g2_direction = "plus"
        g2_pos = max(contig.find(ref2_plus1), contig.find(ref2_plus2))
    elif contig.find(ref2_minus1) != -1 or contig.find(ref2_minus2) != -1 :
        g2_pos = max(contig.find(ref1_minus1), contig.find(ref1_minus2))
    else:
        g2_direction = "error"

    g1_pos_contig = "start"
    g2_pos_contig = "end"
    print(g1_pos, g2_pos, g1_direction, g2_direction)
    print(ref1_plus1)
    print(ref1_plus2)
    print(ref2_minus1)
    print(ref2_minus2)
    if g1_pos > g2_pos :
        g1_pos_contig = "end"
        g2_pos_contig = "start"


    alt1 = ""
    if g1_pos_contig == "start" and g1_direction == "plus" :
        alt1 = "]" + bp2 + "]" + ref2
    elif g1_pos_contig == "start" and g1_direction == "minus" :
        alt1 = "[" + bp1 + "[" + ref2
    elif g1_pos_contig == "end" and g1_direction == "plus" :
        alt1 = ref2 + "[" + bp1 + "["
    elif g1_pos_contig == "end" and g1_direction == "minus" :
        alt1 = ref2 + "]" + bp1 + "]"

    alt2 = ""
    if g2_pos_contig == "start" and g2_direction == "plus" :
        alt2 = "]" + bp2 + "]" + ref1
    elif g2_pos_contig == "start" and g2_direction == "minus" :
        alt2 = "[" + bp2 + "[" + ref1
    elif g2_pos_contig == "end" and g2_direction == "plus" :
        alt2 = ref1 + "[" + bp2 + "["
    elif g2_pos_contig == "end" and g2_direction == "minus" :
        alt2 = ref1 + "]" + bp2 + "]"

    print(alt1)
    print(alt2)

    if g1_pos < g2_pos :
        id = gene1 + "-" + gene2
    else :
        id = gene2 + "-" + gene1

    info = "SVTYPE=Fusion;GENE_NAME=%s;Fusion_support=%s;Gene1_ref_reads=%s;Gene2_ref_reads=%s;Breakpoint1=%s;Breakpoint2=%s" % (
        gene1, fusion_support, gene1_ref_reads, gene2_ref_reads, bp1, bp2
    )
    format = "FS:G1REF:G2REF:GT:GQ"
    data = "%s:%s:%s:./.:." % (fusion_support, gene1_ref_reads, gene2_ref_reads)
    out_line1 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1, start_pos1, id, ref1, alt2, qual, filter, info, format, data)
    fusion_vcf.write(out_line1)
    info = "SVTYPE=Fusion;GENE_NAME=%s;Fusion_support=%s;Gene1_ref_reads=%s;Gene2_ref_reads=%s;Breakpoint1=%s;Breakpoint2=%s" % (
       gene2, fusion_support, gene1_ref_reads, gene2_ref_reads, bp1, bp2
    )
    out_line2 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2, start_pos2, id, ref2, alt1, qual, filter, info, format, data)
    fusion_vcf.write(out_line2)

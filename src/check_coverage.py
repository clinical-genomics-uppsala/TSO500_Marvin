
import subprocess
import sys

#path = "/gluster-storage-volume/data/ref_data/wp1/refFiles_20191010/refFiles/"
path = "/data/ref_data/wp1/refFiles_20191010/refFiles/"
#infiles = ["Mutations_Colon_20171219.csv", "Mutations_Lung_20190211.csv", "Mutations_Melanom_20190211.csv", "Mutations_Ovarial_20170125.csv", "Mutations_Gist_20190211.csv"]
infiles = ["Mutations_Colon_20171219.csv", "Mutations_Lung_20190211.csv", "Mutations_Melanom_20190211.csv", "Mutations_Gist_20190211.csv","Mutations_Ovarial_20170125_header_ok.csv"]

bam_file = sys.argv[1]
vcf_file = open(sys.argv[2])
outfile = open(sys.argv[3], "w")
outfile2 = open(sys.argv[4], "w")


outfile.write("#Chr\tStart_hg19\tEnd_hg19\tGene\tCDS_mut_syntax\tAA_mut_synta\tReport\tcomment\tExon\tAccession_number\tCoverage\tPosition\n")
outfile2.write("#Chr\tStart_hg19\tEnd_hg19\tGene\tCDS_mut_syntax\tAA_mut_synta\tReport\tcomment\tExon\tAccession_number\tCoverage\tPosition")
outfile2.write("\tDP\tRef_DP\tAlt_DP\tAF\tAA_change\tCDS_change\n")


'''Find positions to report and gene regions to analyse'''
inv_pos = {}
gene_regions = []
for infile_name in infiles :
    infile = open(path + infile_name)
    print(infile_name)
    header = True
    prev_gene = ""
    prev_chrom = ""
    first_gene = True
    gene_start_pos = 0
    gene_end_pos = 0
    for line in infile :
        if header:
            header = False
            continue
        lline = line.strip().split("\t")
        Report = lline[6]
        if Report == "region" :
            continue
        chrom = lline[0].split(".")[0].split("0")[-1]
        start_pos = lline[1]
        end_pos = lline[2]
        gene = lline[3]
        pos = int(start_pos)
        while pos <= int(end_pos) :
            inv_pos[chrom + "_" + str(pos)] = lline
            pos += 1
        if first_gene :
            prev_gene = gene
            prev_chrom = chrom
            gene_start_pos = start_pos
            gene_end_pos = end_pos
            first_gene = False
            continue
        if prev_gene != gene :
            gene_regions.append([int(prev_chrom), "chr" + prev_chrom + ":" + gene_start_pos + "-" + gene_end_pos])
            prev_gene = gene
            prev_chrom = chrom
            gene_start_pos = start_pos
            gene_end_pos = end_pos
        else :
            gene_end_pos = end_pos
    gene_regions.append([int(chrom), "chr" + chrom + ":" + gene_start_pos + "-" + gene_end_pos])

'''Remove identical regions'''
gene_regions.sort()
gene_regions_temp = []
for gene_region in gene_regions :
    gene_region = gene_region[1]
    if gene_region not in gene_regions_temp :
        gene_regions_temp.append(gene_region)
gene_regions = gene_regions_temp


'''find all calls in the vcf overlapping hotspots'''
vcf_dict = {}
header = True
for line in vcf_file:
    if header:
        if line[:6] == "#CHROM":
            header = False
        continue
    lline = line.strip().split("\t")
    chrom = lline[0][3:]
    pos = lline[1]
    key = chrom + "_" + pos
    INFO = lline[7].split(";")
    FORMAT = lline[8].split(":")
    DATA = lline[9].split(":")
    if INFO[:3] == "AA=":
        continue
    AD_index = 0
    DP_index = 0
    RD_index = 0
    i = 0
    for f in FORMAT:
        if f == "AD":
            AD_index = i
        if f == "DP":
            DP_index = i
        if f == "RD":
            RD_index = i
        i += 1
    AD = DATA[AD_index].split(",")
    Ref_DP = 0
    Alt_DP = 0
    if len(AD) == 2:
        Ref_DP = AD[0]
        Alt_DP = AD[1]
    else:
        Ref_DP = DATA[RD_index]
        Alt_DP = DATA[AD_index]
    DP = DATA[DP_index]
    AF_index = 0
    i = 0
    for info in INFO:
        if info[:3] == "AF=":
            AF_index = i
        i += 1
    AF = INFO[AF_index][3:]
    VEP = INFO[-1]
    AA_change = VEP.split(":p.")
    if len(AA_change) == 2 :
        AA_change = AA_change[1].split("|")[0]
    else :
        AA_change = ""
    CDS_change = VEP.split(":c.")[1].split("|")[0]
    if key in inv_pos:
        vcf_dict[key] = [DP, Ref_DP, Alt_DP, AF, AA_change, CDS_change]


'''Report all interesting positions (All but region) with coverage < 50'''
depth_dict = {}
for region in gene_regions :
    print(region)
    sample = bam_file.split("/")[-1].split("-ready")[0]
    outfile_name = "DATA/gene_depth_" + sample + ".txt"
    print(sample)
    print(outfile_name)
    subprocess.call("samtools depth -a -r " + region + " " + bam_file + " > " + outfile_name, shell=True)
    depth_file = open(outfile_name)
    for line in depth_file :
        lline = line.strip().split("\t")
        chrom = lline[0][3:]
        pos = lline[1]
        key = chrom + "_" + pos
        if key in inv_pos :
            coverage = int(lline[2])
            if coverage < 50 :
                for info in inv_pos[key] :
                    outfile.write(info + "\t")
                outfile.write(str(coverage) + "\t" + pos + "\n")
            for info in inv_pos[key] :
                outfile2.write(info + "\t")
            outfile2.write(str(coverage) + "\t" + pos)
            if key in vcf_dict:
                for info in vcf_dict[key]:
                    outfile2.write("\t" + str(info))
            outfile2.write("\n")
    depth_file.close()
outfile.close()
outfile2.close()

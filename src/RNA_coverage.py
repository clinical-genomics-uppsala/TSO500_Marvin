
import subprocess
import sys

genes = ["PIK3CA", "MYC", "EML4"]

bam_file = sys.argv[1]
bedfilename = sys.argv[2]
outfile = open(sys.argv[3], "w")
outfile.write("Sample\tGene\tAvg_coverage\t200bp_avg_bins\n")


for gene in genes :
    regions = []
    bedfile = open(bedfilename)
    for line in bedfile :
        if line.find(gene + "_Exon") != -1 and not (line.find("Additional") != -1 or line.find("Fusion") != -1 or line.find("Amp") != -1) :
            lline = line.strip().split("\t")
            regions.append([lline[0], lline[1], lline[2]])
    bedfile.close()
    coverage_sum = 0
    coverage_nr_pos = 0
    coverage_list = []
    for region in regions :
        region = region[0] + ":" + region[1] + "-" + region[2]
        sample = bam_file.split("/")[-1].split(".bam")[0]
        cov_outfile_name = "DATA/RNA_gene_depth_" + sample + ".txt"
        print("samtools depth -a -r " + region + " " + bam_file + " > " + cov_outfile_name)
        subprocess.call("samtools depth -a -r " + region + " " + bam_file + " > " + cov_outfile_name, shell=True)
        depthfile = open(cov_outfile_name)
        for line in depthfile :
            coverage = int(line.strip().split("\t")[2])
            coverage_sum += coverage
            coverage_nr_pos += 1
            coverage_list.append(coverage)
        depthfile.close()
    cov_bins = []
    cov_bin = 0
    i = 0
    for cov in coverage_list :
        if i > 0 and i % 200 == 0 :
            cov_bins.append(cov_bin / 200.0)
            cov_bin = 0
        cov_bin += cov
        i += 1
    if i > 0 and i % 200 == 0 :
        cov_bins.append(cov_bin / 200.0)
    else :
        cov_bins.append(cov_bin / float((i % 200)))
    outfile.write(sample + "\t" + gene + "\t" + str(round(coverage_sum/float(coverage_nr_pos),1)))
    for cov_bin in cov_bins :
        outfile.write("\t" + str(round(cov_bin,1)))
    outfile.write("\n")
outfile.close()

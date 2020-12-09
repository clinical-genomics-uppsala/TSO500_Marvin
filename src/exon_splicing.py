
import subprocess
import sys

bed_file = open(sys.argv[1])
junction_files = sys.argv[2:]

#
#
# bed_file = open("TST500C_manifest.bed")
# junction_files = ["BMS-17SJ.out.tab", "NTRK_13602SJ.out.tab", "R-18-17SJ.out.tab", "R18-295SJ.out.tab", "R-19-135SJ.out.tab", "R19-189SJ.out.tab", "R19-26SJ.out.tab",
#                 "R-19-65SJ.out.tab", "R4_S10SJ.out.tab", "R6_S12SJ.out.tab", "BMS-37SJ.out.tab", "NTRK_29511SJ.out.tab", "R18-226SJ.out.tab", "R-19-120SJ.out.tab",
#                 "R-19-164SJ.out.tab", "R-19-19SJ.out.tab", "R-19-29SJ.out.tab", "R3_S9SJ.out.tab", "R5_S11SJ.out.tab", "R7_S13SJ.out.tab"]

outfile = open("Results/RNA/Exon_skipping/exon_skipping.txt", "w")
outfile.write("Sample\tGene\tChrom\tStart_pos\tEnd_pos\tratio_not_skipped_vs_skipped(<0.05)\tskipped_exons\n")
outfile2 = open("Results/RNA/Exon_skipping/exon_skipping_details.txt", "w")
outfile2.write("Sample\tGene\tChrom\tStart_pos\tEnd_pos\texon\tsplit_reads_in\tsplit_reads_out\tratio_not_skipped_vs_skipped\n")



region_list = []
region_info = []
gene_regions_dict = {}
exon_split_dict = {}
gene_split_dict = {}
i = 0
for line in bed_file :
    lline = line.strip().split("\t")
    chrom = lline[0][3:]
    start_pos = int(lline[1])
    end_pos = int(lline[2])
    region_name = lline[3]
    if region_name.find("MET_Exon") != -1 :
        region_list.append([chrom, start_pos, end_pos, region_name, i])
        gene = region_name.split("_")[0]
        if gene not in gene_regions_dict :
            gene_regions_dict[gene] = []
        gene_regions_dict[gene].append([region_name, str(start_pos), str(end_pos)])
        if chrom not in exon_split_dict :
            exon_split_dict[chrom] = []
        exon_split_dict[chrom].append([start_pos, end_pos, 0, 0, region_name])
        gene = region_name.split("_")[0]
        if gene not in gene_split_dict :
            gene_split_dict[gene] = [0,[],0]
        gene_split_dict[gene][1].append(region_name)
        i += 1


exon_skip_dict = {}
black_list_dict = {}
for jf in junction_files :
    junction_file = open(jf)
    sample = jf.split("/")[-1].split("SJ")[0]
    print(sample)
    '''Exon skipping'''
    for line in junction_file :
        lline = line.strip().split("\t")
        chrom = lline[0]
        start_pos = int(lline[1])
        end_pos = int(lline[2])
        annotated = int(lline[5])
        split_reads = int(lline[6])
        if annotated == 1 and split_reads >= 0 and chrom in exon_split_dict :
            i = 0
            for exon in exon_split_dict[chrom] :
                if start_pos >= exon[0] and start_pos <= exon[1] :
                    exon_split_dict[chrom][i][3] += split_reads
                    gene = exon_split_dict[chrom][i][4].split("_")[0]
                    gene_split_dict[gene][0] += split_reads
                if end_pos >= exon[0] and end_pos <= exon[1] :
                    exon_split_dict[chrom][i][2] += split_reads
                i += 1
    for chrom in exon_split_dict :
        for exon in exon_split_dict[chrom] :
            gene = exon[4].split("_")[0]
            if gene_split_dict[gene][0] > 200 :
                if exon[2] == 0 or exon[3] == 0 :
                    gene_split_dict[gene][2] += 1
    for chrom in exon_split_dict :
        for exon in exon_split_dict[chrom] :
            gene = exon[4].split("_")[0]
            ratio = (exon[2] + exon[3]) / (2 * gene_split_dict[gene][0] / float(len(gene_split_dict[gene][1])))
            outfile2.write(sample + "\t" + gene + "\t" + chrom + "\t" + str(exon[0]) + "\t" + str(exon[1]) + "\t" + exon[4] + "\t" + str(exon[2]) + "\t" + str(exon[3]) + "\t" + str(ratio) + "\n")
            if ratio < 0.05 :
                print(ratio, exon)
                #Sample\tGene\tChrom\tStart_pos\tEnd_pos\tratio_skipped_not_skipped\tskipped_exons\n")
                outfile.write(sample + "\t" + gene + "\t" + chrom + "\t" + str(exon[0]) + "\t" + str(exon[1]) + "\t" + str(ratio) + "\t" + exon[4] + "\n")
    for chrom in exon_split_dict :
        i = 0
        for exon in exon_split_dict[chrom] :
            gene = exon[4].split("_")[0]
            gene_split_dict[gene][0] = 0
            exon_split_dict[chrom][i][2] = 0
            exon_split_dict[chrom][i][3] = 0
            i += 1

out_list = []
for key in exon_skip_dict :
    for sample_exon in exon_skip_dict[key] :
        out_list.append(sample_exon)
out_list.sort()
#for key in exon_skip_dict :
#    for sample_exon in exon_skip_dict[key] :
for sample_exon in out_list :
    outfile.write(sample_exon[0] + "\t" + sample_exon[3][0][3].split("_")[0] + "\t" + str(sample_exon[1]) + "\t" + str(sample_exon[2]) + "\t" + sample_exon[3][0][3])
    for skipped_exon in sample_exon[3][1:] :
        outfile.write("," + skipped_exon[3])
    outfile.write("\n")
    #outfile.write(key.split("__")[0] + "\t" + key.split("__")[1] + "\t" + str(len(exon_skip_dict[key])) + "\n")
outfile.close()
outfile2.close()

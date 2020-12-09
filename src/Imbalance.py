
import subprocess
import sys

bam_files = sys.argv[1:]

gene_files = ["ALK", "NRG1", "RET", "NTRK3", "NTRK1", "NTRK2", "ROS1"]
#gene_files = ["ALK", "NRG1", "RET", "NTRK3"]
outfile = open("Results/RNA/Imbalance/imbalance_all_gene.txt", "w")
outfile2 = open("Results/RNA/Imbalance/imbalance_called_gene.txt", "w")
outfile.write("Sample\tgene\t3'/5'-ratio\tCoverage\tNorm_Coverage\tImbalance\n")
outfile2.write("Sample\tgene\t3'/5'-ratio\tCoverage\tNorm_Coverage\n")

bam_dict = {}
for bam in bam_files :
    print("samtools idxstats " + bam + " | awk -F '\t' '{s+=$3}END{print s}' > nr_reads.txt")
    #subprocess.call("samtools index " + bam, shell=True)
    subprocess.call("samtools idxstats " + bam + " | awk -F '\t' '{s+=$3}END{print s}' > nr_reads.txt", shell=True)
    nr_bam_file = open("nr_reads.txt")
    nr_reads = 0
    for line in nr_bam_file :
        nr_reads = int(line.strip()) / 10000000.0
    bam_dict[bam] = nr_reads
    nr_bam_file.close()

subprocess.call("mkdir exon_coverage", shell = True)
exon_result_dict = {}
for bam in bam_files :
    exon_result_dict[bam] = {}
    print(bam)
    print(bam_dict[bam])
    #sample = bam.split("/")[-1].split("Aligned")[0]
    sample = bam.split("/")[-1].split(".bam")[0]
    for gene in gene_files :
        exon_result_dict[bam][gene] = []
        print(gene)
        gene_file = open("DATA/" + gene + "_exons.txt")
        exon_list = []
        for line in gene_file :
            lline = line.strip().split("\t")
            chrom = lline[0]
            start_pos = lline[1]
            end_pos = lline[2]
            exon = lline[3]
            exon_list.append([exon, chrom, start_pos, end_pos])
            #print(exon, chrom, start_pos, end_pos)
            print("samtools depth -aa -r " + chrom + ":" + str(start_pos) + "-" + str(end_pos) + " " + bam + " > exon_coverage/" + exon + ".txt")
            subprocess.call("samtools depth -aa -r " + chrom + ":" + str(start_pos) + "-" + str(end_pos) + " " + bam + " > exon_coverage/" + exon + ".txt", shell=True)
        exon_result_dict[bam][gene].append(exon_list)
        avg_depth_list = []
        i = 0
        for exon in exon_list :
            exon_file = open("exon_coverage/" + exon[0] + ".txt")
            avg_depth = 0
            nr_pos = 0
            for line in exon_file :
                lline = line.strip().split("\t")
                chrom = lline[0]
                pos = lline[1]
                avg_depth += int(lline[2])
                nr_pos += 1
            avg_depth /= float(nr_pos)
            avg_depth_list.append(avg_depth)
            exon_result_dict[bam][gene][0][i].append(avg_depth)
            i += 1
            exon_file.close()

        outfile.write(sample + "\t" + gene + "\t")
        gene_start = 1
        gene_end = 1
        coverage = 0
        i = 0
        '''ALK'''
        if gene == "ALK" :
            for d in avg_depth_list :
                if i <= 2 :
                    gene_start += d / bam_dict[bam]
                elif i > 18 :
                    gene_end += d / bam_dict[bam]
                if i > 18 :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 10.0
            gene_start /= 3.0
            coverage /= 10.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            #if imbalance > 3.0 and coverage > 40.0 :
            if imbalance > 25.0 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''NRG1'''
        if gene == "NRG1" :
            for d in avg_depth_list :
                if i <= 2 :
                    gene_start += d / bam_dict[bam]
                elif i >= 5 :
                    gene_end += d / bam_dict[bam]
                if i >= 5 :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 8.0
            gene_start /= 3.0
            coverage /= 8.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            if imbalance > 3.0 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''RET'''
        if gene == "RET" :
            for d in avg_depth_list :
                if i <= 2 :
                    gene_start += d / bam_dict[bam]
                elif i >= 11 :
                    gene_end += d / bam_dict[bam]
                if i >= 11 :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 9.0
            gene_start /= 3.0
            coverage /= 9.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            if imbalance > 3.0 and coverage > 20.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''NTRK3'''
        if gene == "NTRK3" :
            for d in avg_depth_list :
                if i <= 2 :
                    gene_start += d / bam_dict[bam]
                elif i >= 14 :
                    gene_end += d / bam_dict[bam]
                if i >= 14  :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 6.0
            gene_start /= 3.0
            coverage /= 6.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            #if imbalance > 3.0 and coverage > 40.0 :
            if imbalance > 1.0 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''NTRK1'''
        if gene == "NTRK1" :
            for d in avg_depth_list :
                if i >= 11 :
                    gene_end += d / bam_dict[bam]
                elif i <= 2 :
                    gene_start += d / bam_dict[bam]
                if i >= 11  :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 6.0
            gene_start /= 3.0
            coverage /= 6.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            if imbalance > 3.0 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''NTRK2'''
        if gene == "NTRK2" :
            for d in avg_depth_list :
                if i >= 14 :
                    gene_end += d / bam_dict[bam]
                elif i <= 2 :
                    gene_start += d / bam_dict[bam]
                if i >= 14  :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 7.0
            gene_start /= 2.0
            coverage /= 7.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            #if imbalance > 3.0 and coverage > 50.0 :
            if imbalance > 1.0 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
        '''ROS1'''
        if gene == "ROS1" :
            for d in avg_depth_list :
                if i <= 2 :
                    gene_start += d / bam_dict[bam]
                elif i >= 34 :
                    gene_end += d / bam_dict[bam]
                if i >= 34  :
                    coverage += d / bam_dict[bam]
                i += 1
            gene_end /= 9.0
            gene_start /= 2.0
            coverage /= 9.0
            imbalance = gene_end / gene_start
            outfile.write(str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\t")
            if imbalance > 1.1 and coverage > 40.0 :
                outfile.write("Imbalance" + "\n")
                outfile2.write(sample + "\t" + gene + "\t" + str(round(imbalance,1)) + "\t" + str(round(coverage*bam_dict[bam],1)) + "\t" + str(round(coverage,1)) + "\n")
            else :
                outfile.write("" + "\n")
outfile.close()
outfile2.close()

for gene in gene_files :
    outfile3 = open("Results/RNA/Imbalance/imbalance_" + gene + "_exons.txt", "w")
    outfile3.write("Sample")
    for exon in exon_result_dict[bam][gene][0] :
        outfile3.write("\t" + exon[0])
    outfile3.write("\n")
    for bam in bam_files :
        sample = bam.split("/")[-1].split("Aligned")[0]
        i = 0
        for exon in exon_result_dict[bam][gene][0] :
            if i == 0 :
                outfile3.write(sample)
            outfile3.write("\t" + str(round(exon[4],1)))
            i += 1
        outfile3.write("\n")
    outfile3.close()


import sys

outfile = open(sys.argv[1], "w")
infiles = sys.argv[2:]

#outfile = open("Bam_stats.txt", "w")
#infiles = ["R20-16/QC/RSeQC_bam_stat.txt", "R20-15/QC/RSeQC_bam_stat.txt"]

first = True
bam_stat = {}
columns = []
for infile_name in infiles :
    sample = infile_name.split("/")[2]
    bam_stat[sample] = []
    infile = open(infile_name)
    header = True
    for line in infile :
        if header :
            if line.find("Total") != -1 :
                header = False
            else :
                continue
        column_name = ""
        count = 0
        if line.strip() == "" :
            continue
        for s in line.strip() :
            if s == ":" :
                break
            column_name += s
        count = int(line.strip().split(" ")[-1].split(":")[-1])
        if first :
            columns.append(column_name)
        bam_stat[sample].append(str(count))
    first = False
    infile.close()

outfile.write("Measure")
for sample in bam_stat :
    outfile.write("\t" + sample)
outfile.write("\n")
i = 0
for column in columns :
    outfile.write(column)
    for sample in bam_stat :
        outfile.write("\t" + bam_stat[sample][i])
    outfile.write("\n")
    i += 1
outfile.close()

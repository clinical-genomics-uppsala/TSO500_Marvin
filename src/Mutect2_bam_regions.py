
import gzip
import subprocess
import sys

vcf_filename = sys.argv[1]
bam_filename = sys.argv[2]
sample = bam_filename.split("/")[-1].split("-ready")[0]
flank_size = 1000
chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]

with gzip.open(vcf_filename, "rt") as infile:
    header = True
    pos_list = []
    for line in infile :
        if header :
            if line.find("#CHROM") != -1 :
                header = False
            continue
        lline = line.split("\t")
        ref = lline[3]
        alt = lline[4]
        if len(ref) + len(alt) <= 2 :
            continue
        chrom = lline[0]
        pos = int(lline[1])
        #print(region[0] + "\t" + str(region[1]) + "\t" + str(region[2]))
        if pos_list != [] :
            prev_chrom = pos_list[-1][0]
            prev_start_pos = pos_list[-1][1]
            prev_stop_pos = pos_list[-1][2]
            if prev_chrom == chrom and pos - flank_size <= prev_stop_pos :
                pos_list[-1][2] = pos + flank_size
            else :
                pos_list.append([chrom, pos - flank_size, pos + flank_size])
        else :
            pos_list.append([chrom, pos - flank_size, pos + flank_size])

for chrom in chrom_list :
    bed_file = "Mutect2/" + sample + "." + str(chrom) + ".bed"
    outbedfile = open(bed_file, "w")
    for region in pos_list :
        if region[0][3:] == str(chrom) :
            outbedfile.write(region[0] + "\t" + str(region[1]) + "\t" + str(region[2]) + "\n")
    outbedfile.close()
    #command_line = "samtools view " + bam_filename + " -b -L " + bed_file + " > Mutect2/" + sample + "-ready." + str(chrom) + ".indel.bam"
    #command_line += " && samtools index Mutect2/" + sample + "-ready." + str(chrom) + ".indel.bam"
    #print(command_line)
    #subprocess.call(command_line, shell=True)


import sys

refseq = open(sys.argv[1])
infile = open(sys.argv[2])

gene_list = []
gene_dict = {}
header = True
header_line = ""
i = 0
for line in infile :
    if header :
        header = False
        header_line = "#fusion\t" + line[1:]
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    chr1 = lline[4].split(":")[0]
    chr2 = lline[5].split(":")[0]
    pos1 = int(lline[4].split(":")[1])
    pos2 = int(lline[5].split(":")[1])
    gene_list.append([[gene1, chr1, pos1], [gene2, chr2, pos2], line, "", "", "", ""])
    if gene1 not in gene_dict:
        gene_dict[gene1] = [[i,3]]
    else :
        gene_dict[gene1].append([i,3])
    if gene2 not in gene_dict:
        gene_dict[gene2] = [[i,4]]
    else :
        gene_dict[gene2].append([i,4])
    i += 1
    #print(gene_list)
infile.close()

header = True
for line in refseq :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene_name = lline[9]
    if gene_name in gene_dict :
        i_list = gene_dict[gene_name]
        for i in i_list :
            if gene_list[i[0]][i[1]] == "" :
                transcript = lline[0]
                exonstarts = lline[7].split(",")
                exonends = lline[8].split(",")
                exon_nr = 0
                pos = gene_list[i[0]][i[1]-3][2]
                for startpos in exonstarts[:-1] :
                    print(pos,startpos,exonends[exon_nr])
                    if pos >= int(startpos) and pos <= int(exonends[exon_nr]) :
                        gene_list[i[0]][i[1]] = str(exon_nr)
                        gene_list[i[0]][i[1]+2] = transcript
                        break
                    exon_nr += 1

outfile = open(sys.argv[3], "w")
outfile.write(header_line)
for fusion_info in gene_list :
    #print(fusion_info)
    fusion = fusion_info[0][0] + "_" + fusion_info[5] + "_exon" + fusion_info[3] + ":" + fusion_info[1][0] + "_" + fusion_info[6] + "_exon" + fusion_info[4]
    outfile.write(fusion + "\t" + fusion_info[2])
outfile.close()


localrules: SampleSheet_TSO500

rule SampleSheet_TSO500:
    input:
        Sample_sheet = config["Sample_sheet"]
    output:
        Sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    params:
        DNA_samples = expand(config["DNA_Samples"])
    run:
        SS = open(output.Sample_sheet, "w")
        infile = open(input.Sample_sheet)
        header = True
        for line in infile :
            if header :
                if line[:5] == "Lane," :
                    header = False
                    SS.write("Lane,Sample_ID,index,index2,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,I5_Index_ID,Project,Description\n")
                else :
                    SS.write(line)
                continue
            sample = line.split(",")[1]
            if sample in params.DNA_samples :
                lline = line.strip().split(",")
                SS.write(lline[0])
                for column in lline[1:-1] :
                    if column == "DNA" :
                        column = ""
                    SS.write("," + column)
                SS.write("\n")
        SS.close()
        infile.close()

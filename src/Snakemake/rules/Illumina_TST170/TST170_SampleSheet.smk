
localrules: SampleSheet_TST170


rule SampleSheet_TST170:
    input:
        Sample_sheet = config["Sample_sheet"]
    output:
        Sample_sheet = "SampleSheet.csv"
    params:
        RNA_samples = expand(config["RNA_Samples"]),
        Runfolder = config["Runfolder"]
    run:
        import subprocess
        SS = open(output.Sample_sheet, "w")
        infile = open(input.Sample_sheet)
        state = 0
        for line in infile :
            if state == 0 :
                if line[:7] == "[Reads]" :
                    state = 1
                    SS.write("[Manifests],,,,,,,,\n")
                    SS.write("PoolDNA,MixDNA_Manifest.txt,,,,,,,\n")
                    SS.write("PoolRNA,MixRNA_Manifest.txt,,,,,,,\n")
                SS.write(line)
                continue
            elif state == 1 :
                if line[:5] == "Lane," :
                    SS.write("Lane,Sample_ID,index,index2,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,I5_Index_ID,Project,Description,Manifest\n")
                    state = 2
                    continue
                else :
                    SS.write(line)
            elif state == 2 :
                sample = line.split(",")[1]
                if sample in params.RNA_samples :
                    lline = line.strip().split(",")
                    SS.write(lline[0])
                    for column in lline[1:-1] :
                        SS.write("," + column)
                    SS.write(",PoolRNA\n")
        SS.close()
        infile.close()
        subprocess.call("cp " + output.Sample_sheet + " " + params.Runfolder, shell=True)

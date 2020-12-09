
localrules: all, Create_TSO500_yaml

rule all:
    input:
        TSO500_yaml = "TSO500.yaml",
        TC = "DATA/Pathological_purity_BMS_validation.txt"


rule Create_TSO500_yaml:
    input:
        run_info = "RunParameters.xml",
        config = "Config/Pipeline/configdefaults201012.yaml"
    output:
        TSO500_yaml = "TSO500.yaml",
        TC = "DATA/Pathological_purity_BMS_validation.txt"
    run:
        import glob
        import os
        import subprocess
        subprocess.call("cp " + input.config + " " + output.TSO500_yaml, shell=True)
        run_folder_name = ""
        run_info_file = open(input.run_info)
        for line in run_info_file :
            if line.find("<RunID>") != -1 :
                run_folder_name = line.split(">")[1].split("<")[0]
        state = 0
        DNA_sample_list = []
        RNA_sample_list = []
        #KG_runname = os.getcwd().split("/")[-1]
        i = 1
        sample_sheet_name = glob.glob("*heet.csv")
        if len(sample_sheet_name) > 1 :
            print("Error: Something wrong with the sample sheet name!")
            quit()
        sample_sheet_name = sample_sheet_name[0]
        infile = open(sample_sheet_name)
        for line in infile:
            if state == 0 :
                if line.find("Experiment Name") != -1 :
                    state = 1
            elif state == 1 :
                if line.find("Lane,Sample_ID") != -1 :
                    state = 2
            elif state == 2 :
                lline = line.strip().split(",")
                sample = lline[1]
                if sample.find(" ") != -1 :
                    print("incorrect sample name: " + sample)
                    quit()
                sample_type = lline[9]
                if sample_type == "DNA" :
                    TC = lline[11]
                    DNA_sample_list.append([sample, i, TC])
                elif sample_type == "RNA" :
                    RNA_sample_list.append([sample, i])
                else :
                    print("Error: wrong sample type: " + sample_type)
                    quit()
                i += 1
        outfile = open(output.TSO500_yaml, "a")
        outfile2 = open(output.TC, "w")
        #outfile.write("Runfolder: /projects/wp1/nobackup/ngs/klinik/INBOX/" + KG_runname + "/\n\n")
        outfile.write("Runfolder: /projects/wp1/nobackup/ngs/klinik/INBOX/" + run_folder_name + "/\n\n")
        #outfile.write("Runfolder: " + run_folder_name + "/\n\n")
        #outfile.write("Outfolder: /projects/wp1/nobackup/ngs/klinik/OUTBOX/" + KG_runname + "/\n\n")
        outfile.write("Sample_sheet: " + sample_sheet_name + "\n\n")
        if len(DNA_sample_list) == 0 :
            outfile.write("DNA_Samples: No DNA\n")
        else :
            outfile.write("DNA_Samples:\n")
            for sample in DNA_sample_list :
                outfile.write("  " + sample[0] + ": \"S" + str(sample[1]) + "\"\n")
                outfile2.write(sample[0] + "-ready\t" + sample[2] + "\n")
        if len(RNA_sample_list) == 0 :
            outfile.write("\nRNA_Samples: No RNA\n")
        else :
            outfile.write("\nRNA_Samples:\n")
            for sample in RNA_sample_list :
                outfile.write("  " + sample[0] + ": \"S" + str(sample[1]) + "\"\n")
        outfile.close()
        outfile2.close()

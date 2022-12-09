# Howto run TSO500 manually 
## 1 - Copy data
### 1.1 Copy entire runfolder from hospital server to Marvin
```
rsync -Prv {runfolder} user@server:/data/INBOX/
```
### 1.2 - Clone git and copy needed runfiles to analysis folder
runfolder is name from sequencing machine. KP_runname is the name given for this run. 
```
cd /data/{runfolder}/
mkdir /scratch/{KP_runname}/
git clone https://github.com/clinical-genomics-uppsala/TSO500_Marvin.git /scratch/{KP_runname}/
rsync -Prv *sheet.csv *.xml /scratch/{KP_runname}/
```

## 2 Create TSO500.yaml
### Snakemake
```
screen -S TSO500_{DATE}
cd /scratch/{KP_runname}/
module add slurm-drmaa
module add singularity
module add snakemake
snakemake -p -j 1 --drmaa "-A wp1 -p core -n 1 -t 2:00:00 "  -s ./src/Snakemake/rules/TSO500_yaml/TSO500_yaml.smk
```
 
# Troubleshooting
## ERROR in Samplesheet creation?
Doublecheck samplesheet:
 - No blank spaces allowed
 - TC column is required
 - DNA or RNA must be filled in in the Project column
 - First column should be empty

When corrected rerun previous snakemake command
```
snakemake -p -j 80 --drmaa "-A wp1 -p core -n {cluster.n} -t {cluster.time} -N 1-1"  -s ./TSO500.smk --restart-times 2 --use-singularity --singularity-prefix /data/singularity_cache/ --singularity-args "-e --bind /data --bind /projects --bind /scratch --cleanenv " --cluster-config Config/Slurm/cluster.json
```

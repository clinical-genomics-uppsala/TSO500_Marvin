

rule check_DNA_coverage:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.AD20.vep.multiplebp.vcf",
    output:
        coverage = "Results/DNA/{sample}/QC/Low_coverage_positions.txt",
        coverage2 = "Results/DNA/{sample}/QC/All_coverage_positions.txt",
    run:
        import subprocess
        subprocess.call("python src/check_coverage.py " + input.bam + " " + input.vcf + " " + output.coverage + " " + output.coverage2, shell=True)

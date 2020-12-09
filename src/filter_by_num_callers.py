import gzip
import io
import re
import argparse


parser = argparse.ArgumentParser(description = "This script filter vcf file using a file with defined variants.")
parser.add_argument('-v', '--variantfile', help = 'VCF-file', type = str, required = True)
parser.add_argument('-d', "--deletions", help ="Include deletions with one caller support", action="count", default=0)
#parser.add_argument('-o', '--blacklist', help = 'Blacklist', type = str, required = False)


args = parser.parse_args()

with gzip.open(args.variantfile,'r') as lines:
	with io.TextIOWrapper(lines, encoding='utf-8') as enc:
		for line in enc:
			if line.startswith('#'):
				print(line, end='')	
			else:
				data = re.search('^chr[XYMT0-9]+\t\d+\t[.|sr0-9;]+\t([ACGT]+)\t([ACGT,]+)\t[0-9.]+\t\w+\t.*CALLERS=([a-z,2]+);.+', line, re.IGNORECASE)
				ref = data.group(1)
				var = data.group(2)
				callers =  data.group(3).split(",")
				if len(callers) > 1:
                                        print(line, end='')
				elif args.deletions and (len(var) > 1 or len(ref) > 1) and ("mutect2" in callers or "vardict" in callers):
					print(line, end='')


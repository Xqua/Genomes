#!/usr/bin/env python 

from Bio import AlignIO
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--stockholm_file", dest="stockholmpath", default="None",
                  help="[Required] Location of the alignement file")
parser.add_option("-o", "--output", dest="fastapath", default="alignment.fasta",
                  help="Location of the fasta output alignement file")

(options, args) = parser.parse_args()

fastapath = options.fastapath
stockholmpath = options.stockholmpath
if stockholmpath == "None":
	print "File must be provided"
	sys.exit(1)


# f = open(stockholmpath)
# lines = f.readlines()
# names = {}

# for line in lines:
# 	if "#=GS" in line:
# 		ID = line.split(' ')[1]
# 		genome = line.split('|')[-1].strip()
# 		names[ID] = genome
# f.close()



out = open(fastapath,'w')
align = AlignIO.read(stockholmpath, "stockholm")
seqs = list(align)

for seq in seqs:
	title = seq.annotations['accession']
	title += '|'
	title += seq.description.replace('[subseq from]','')
	title += '|%s-%s' % (seq.annotations['start'], seq.annotations['end'])
	title = title.replace(' ','')
	out.write('>' + title + '\n')
	out.write(str(seq.seq) + '\n')
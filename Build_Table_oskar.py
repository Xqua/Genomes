#!/usr/bin/env python

from Bio import AlignIO
from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--stockholm_file1", dest="stockholmpath_osk", default="None",
                  help="[Required] Location of the oskar alignement file")
parser.add_option("-j", "--stockholm_file2", dest="stockholmpath_lotus", default="None",
                  help="[Required] Location of the lotus alignement file")
parser.add_option("-k", "--stockholm_file3", dest="stockholmpath_sgnh", default="None",
                  help="[Required] Location of the sgnh alignement file")
parser.add_option("-d", "--nrdb_file", dest="nrdbpath", default="None",
                  help="[Required] Location of the database used to produce the alignement")
parser.add_option("-o", "--output", dest="fastapath", default="hmmer_result.fasta",
                  help="Location of the fasta output file")
parser.add_option("-t", "--TSA",
                  action="store_true", dest="tsadb", default=False,
                  help="IS the file from the TSA DB ?")


(options, args) = parser.parse_args()

fastapath = options.fastapath
stockholmpath_osk = options.stockholmpath_osk
stockholmpath_lotus = options.stockholmpath_lotus
stockholmpath_sgnh = options.stockholmpath_sgnh
nrdbpath = options.nrdbpath
tsadb = options.tsadb

if stockholmpath_osk == "None" or stockholmpath_lotus == "None" or  stockholmpath_sgnh == "None" or  nrdbpath == "None":
	print "Stockholm and NRDB File must all be provided"
	sys.exit(1)


out = open(fastapath,'w')
align_osk = AlignIO.read(stockholmpath_osk, "stockholm")
align_lotus = AlignIO.read(stockholmpath_lotus, "stockholm")
align_sgnh = AlignIO.read(stockholmpath_sgnh, "stockholm")
seqs_osk = list(align_osk)
seqs_lotus = list(align_lotus)
seqs_sgnh = list(align_sgnh)
nrdbhandle = SeqIO.parse(nrdbpath, 'fasta') 

print "Creating genome list ... wait..."
all_genome = {}
for seq in nrdbhandle:
	desc = seq.description
	genome = desc.split('|')[-1].strip()
	geneID = desc.split('|')[-2].strip()
	if not all_genome.has_key(genome):
		all_genome[genome] = [0,0,0,0]
print "Done"

def find_match(seqs):
	match = {}
	for seq in seqs:
		title = seq.annotations['accession']
		title += '|'
		title += seq.description.replace('[subseq from]','')
		title += '|%s-%s' % (seq.annotations['start'], seq.annotations['end'])
		title = title.replace(' ','')
		snap = False
		augu = False
		s = title.split('|')
		print s
		if tsadb:
			geneID = s[4]
			genome = geneID[0:4]
			frame = s[-2]
		else:
			for a in s:
				if 'snap' in a:
					snap = True
				else:
					augu = True
			if snap:
				geneID = s[2]
				genome = s[4]
			elif augu:
				geneID = s[3]
				genome = s[4]
			else:
				print "WTF ?"
				sys.exit(1)
		if not match.has_key(genome):
			match[genome] = []
		if tsadb:
			if [geneID, frame] not in match[genome]:
				match[genome].append([geneID, frame])
		else:
			if geneID not in match[genome]:
				match[genome].append(geneID)
		# print "Added genome %s gene %s to the retrieval list" % (genome, geneID)


print "extracting matches from stockholm files"
match_osk = find_match(seqs_osk)
match_lotus = find_match(seqs_lotus)
match_sgnh = find_match(seqs_sgnh)

for genome in match_osk.keys():
	all_genome[genome][0] = len(match_osk[genome])

for genome in match_lotus.keys():
	all_genome[genome][1] = len(match_lotus[genome])

for genome in match_sgnh.keys():
	all_genome[genome][2] = len(match_sgnh[genome])

for genome1 in match_sgnh.keys():
	for genome2 in match_lotus.keys():
		if genome1 == genome2:
			for gene1 in match_sgnh[genome1]:
				for gene2 in match_lotus[genome2]:
					if tsadb:
						if gene1[0] == gene2[0] and gene1[1] == gene2[1]:
							all_genome[genome1][3] += 1
					else:
						if gene1 == gene2:
								all_genome[genome1][3] += 1
print "Done"
print "Writing results"

res = []
for genome in all_genome:
	res.append([genome, str(all_genome[genome][0]), str(all_genome[genome][1]), str(all_genome[genome][2]), str(all_genome[genome][3])])

out.write('GenomeID\tOskar_Hits\tLOTUS_Hits\tSGNH_Hits\tLOTUS_and_SGNH\n')
for i in res:
	t = '\t'.join(i)
	out.write(t + '\n')
out.close()
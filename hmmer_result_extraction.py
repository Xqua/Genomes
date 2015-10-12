#!/usr/bin/env python

from Bio import AlignIO
from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--stockholm_file", dest="stockholmpath", default="None",
                  help="[Required] Location of the alignement file")
parser.add_option("-d", "--nrdb_file", dest="nrdbpath", default="None",
                  help="[Required] Location of the database used to produce the alignement")
parser.add_option("-o", "--output", dest="fastapath", default="hmmer_result.fasta",
                  help="Location of the fasta output file")
parser.add_option("-t", "--TSA",
                  action="store_true", dest="tsadb", default=False,
                  help="IS the file from the TSA DB ?")


(options, args) = parser.parse_args()

fastapath = options.fastapath
stockholmpath = options.stockholmpath
nrdbpath = options.nrdbpath
tsadb = options.tsadb
if stockholmpath == "None" or nrdbpath == "None":
	print "Stockholm and NRDB File must be provided"
	sys.exit(1)


out = open(fastapath,'w')
align = AlignIO.read(stockholmpath, "stockholm")
seqs = list(align)
nrdbhandle = SeqIO.parse(nrdbpath, 'fasta') 

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
	print "Added genome %s gene %s to the retrieval list" % (genome, geneID)

towrite = []
for seq in nrdbhandle:
	desc = seq.description
	if tsadb:
		geneID = desc.split('|')[4].strip()
		genome = geneID[0:4]
		frame = desc.split('|')[-1].strip()
		# print geneID, genome, frame
		if genome in match.keys():
			# print match[genome]
			for g in match[genome]:
				if geneID in g[0] and frame in g[1]:
					print "Found genome %s gene %s frame %s, saving it" % (genome, geneID, frame)
					towrite.append(seq)
	else:
		genome = desc.split('|')[-1].strip()
		geneID = desc.split('|')[-2].strip()
		if genome in match.keys():
			if geneID in match[genome]:
				print "Found genome %s gene %s, saving it" % (genome, geneID)
				towrite.append(seq)

SeqIO.write(towrite, out, 'fasta')
out.close()
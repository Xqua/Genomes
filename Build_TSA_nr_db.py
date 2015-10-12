#!/usr/bin/env python
# This scripts has been built by Leo Blondel on 9/14/2015
# This script reads from GCA accessed NCBI genomes (see genome_sucking.py)
# Takes all the fna files (DNA sequences) along with the gff (Annotations)
# and builds a Fasta file containing the coresponding protein sequences.

import os, sys
from optparse import OptionParser
import gzip
from Bio import SeqIO

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-o", "--output", dest="dbPath", default="hmmer_db.fasta",
                  help="Location of the output database file")
parser.add_option("-6", "--6frames",
                  action="store_true", dest="doFrames", default=False,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
if genomeFolder == "None":
	print "Genome Folder must be provided"
	sys.exit(1)

doFrames = options.doFrames
dbPath = options.dbPath

gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


def sixFrameTranslate(DNA):
	res = []
	sequences = []
	for frameshift in [0,1,2]:
		sequences.append(DNA[frameshift:])
		sequences.append(ReverseComplement(DNA)[frameshift:])
	for seq in sequences:
		prot = Translate(seq)
		res.append(prot)
	return res

def Translate(seq):
	seq = seq.upper()
	prot = ""
	for i in range(len(seq)/3):
		codon = seq[i*3:i*3+3]
		if not "N" in codon:
			prot += gencode[codon]
	return prot

def ReverseComplement(DNA):
	RCDNA = ""
	for letter in DNA:
		if letter == "A":
			RCDNA += "T"
		elif letter == "T":
			RCDNA += "A"
		elif letter == "C":
			RCDNA += "G"
		elif letter == "G":
			RCDNA += "C"
		elif letter == "N":
			RCDNA += "N"
	return RCDNA

def readFASTA(f):
	"takes FASTA file handle and reads it into a dictionary"
	data = {}
	for line in f.readlines():
		l = line.strip()
		if l[0] == ">":
			t = l[1:]
			data[t] = ""
		else:
			data[t] += l
	return data

out = open(dbPath,'w')
list_genomes = os.listdir(genomeFolder)
for genome in list_genomes:
	print "Opening Genome: %s" % genome
	path = os.path.join(genomeFolder, genome, genome+"_genomic.fna.gz")
	print path
	g = gzip.open(path)
	fasta = readFASTA(g)
	g.close()
	for ID in fasta.keys():
		if doFrames:
			prots = sixFrameTranslate(fasta[ID])
			for prot in prots:
				out.write('>%s|%s\n'%(genome, ID))
				out.write(prot + '\n')
		else:
			prot = Translate(fasta[ID])
			out.write('>%s|%s\n'%(genome, ID))
			out.write(prot + '\n')
	print "genome %s: Added !" % genome


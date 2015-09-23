#!/usr/bin/env python

from subprocess import Popen
import subprocess
from optparse import OptionParser
import sys, os

def get_protein_list(path):
	f = open(path)
	lines = f.readlines()
	f.close()
	proteins = {}
	rec  = False
	for line in lines:
		if "# start gene" in line:
			GID = line.split(' ')[3].strip()
		if "# end gene" in line:
			rec = False
			proteins[GID] = sequence
		if "# protein sequence" in line:
			if ']' in line:
				sequence = line.split('[')[1].split(']')[0]
			else:
				rec = True
				sequence = line.split('[')[1].strip()
		if rec:
			if ']' in line:
				sequence += line.split(' ')[1].split(']')[0]
			else:
				sequence += line.split(' ')[1].strip()
	return proteins

def prot2fasta(proteins, path):
	f = open(path,'w')
	GIDs = proteins.keys()
	GIDs.sort()
	for GID in GIDs:
		f.write('> '+GID+'\n')
		f.write(proteins[GID] + '\n')
	f.close()

def Merge(file1, file2, outPath):
	s = open(file1)
	a = open(file2)
	f = open(outPath, 'w')
	f.write(s.read())
	f.write('\n')
	f.write(a.read())
	f.close()

def uClust(inputPath, ucpath, outputPath, precision):
	cmd = ['uclust',
			 '--sort',
			 inputPath,
			 '--output',
			 os.path.join(ucpath, 'seqs_sorted.fasta')]
	P = Popen(cmd)
	ret = P.wait()

	cmd = ['uclust',
			 '--input',
			 os.path.join(ucpath, 'seqs_sorted.fasta'),
			 '--uc',
			 os.path.join(ucpath,'results.uc'),
			 '--id',
			 precision]
	P = Popen(cmd)
	ret = P.wait()

	cmd = ['uclust',
	 	  '--input',
		  inputPath,
		  '--uc2fasta',
		  os.path.join(ucpath,'results.uc'),
		  '--types',
		  'S',
		  '--output',
		  outputPath]
	P = Popen(cmd)
	ret = P.wait()

def Annotate(nrPath, genome):
	i = open(nrPath)
	lines = i.readlines()
	i.close()
	o = open(nrPath,'w')
	for line in lines:
		if ">" in line:
			o.write(line.strip() + ' | ' + genome + '\n')
		else:
			o.write(line)
	o.close()

def WriteToDB(nrdbpath, dbpath):
	db = open(dbpath,'a')
	f = open(nrdbpath)
	db.write(f.read() + '\n')
	db.close()



parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-d", "--database_path", dest="db_path", default="./nrdb.fasta",
                  help="Location of the output file that will contain the non redundant database")
parser.add_option("-a", "--append",
                  action="store_true", dest="app", default=False,
                  help="Append to the database (This will not erase the DB)")
parser.add_option("-p", "--precision", dest="uclust_precision", default="0.90",
                  help="Overlap Threshold for uClust to call two sequence identical (uClust will keep the longer sequence of the two)")


(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
if genomeFolder == "None":
	print "Genome Folder must be provided.\n -h for help"
	sys.exit(1)
db_path = options.db_path
app = options.app
if not app:
	db = open(db_path,'w')
	db.close()

threshold = options.uclust_precision

list_genomes = os.listdir(genomeFolder)
for genome in list_genomes:
	path = os.path.join(genomeFolder, genome, genome+"_genomic.fna")
	if os.path.isfile(path):
		if os.path.isfile(path+".SNAPoutput") and os.path.isfile(path+'_augustus.gff'):
			print "Opening Genome: %s" % genome
			pr = get_protein_list(path+'_augustus.gff')
			prot2fasta(pr, path+'_augustus.fasta')
			Merge(path+".SNAPoutput", path+'_augustus.fasta', path+'.merged.fasta')
			uClust(path+'.merged.fasta',os.path.join(genomeFolder, genome),  path+'.nrdb.fasta', threshold)
			Annotate(path+'.nrdb.fasta', genome)
			WriteToDB(path+'.nrdb.fasta', db_path)
			print "genome %s: Added to non redundant database !" % genome
		else:
			print "Genome %s does not contain augustus and snap file" % genome
	else:
		print "no Genome found ! O.o Weirrrrddddd"	
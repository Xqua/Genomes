#!/usr/bin/env python
# This scripts has been built by Leo Blondel on 9/14/2015
# This script reads from GCA accessed NCBI genomes (see genome_sucking.py)
# Takes all the fna files (DNA sequences) along with the gff (Annotations)
# and builds a Fasta file containing the coresponding protein sequences.

import os, sys
from optparse import OptionParser
from subprocess import Popen
import subprocess
import time

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-c", "--cpu", dest="cpu", default=1,
                  help="Number of CPU to use for the analysis. Each core is used to launch one instance of SNAP or Augustus.")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
if genomeFolder == "None":
	print "Genome Folder must be provided"
	sys.exit(1)

try:
	nb_CPU = int(options.cpu)
except:
	print "Wrong number of CPUs !"
	sys.exit(1)

cmdToRun = []
list_genomes = os.listdir(genomeFolder)
for genome in list_genomes:
	print "Opening Genome: %s" % genome
	path = os.path.join(genomeFolder, genome, genome+"_genomic.fna.gz")
	print path
	if os.path.isfile(path):
		if not os.path.isfile(os.path.join(genomeFolder, genome, genome+"_genomic.fna")):
			P = Popen(['gunzip',path])
			ret = P.wait()
			if ret != 0:
				print "Error Gunzipping !"
	path = os.path.join(genomeFolder, genome, genome+"_genomic.fna")
	if os.path.isfile(path):
		# geneMarkCMD = ["gmhmme3", "-m", "/home/lblondel/Software/genemark_hmm_euk_linux_64/ehmm/d_melanogaster.mod", "-p", "-n",path]
		# cmdToRun.append(geneMarkCMD)
		if not os.path.isfile(path+".SNAPoutput"):
			SNAPCMD = ['snap','fly',path,'-aa',path+".SNAPoutput",'-gff']
			cmdToRun.append(SNAPCMD)
		if not os.path.isfile(path+'_augustus.gff'):
			AugustusCMD = ["augustus",'--gff3=on','--progress=true','--outfile='+path+'_augustus.gff','--species=fly',path]
			cmdToRun.append(AugustusCMD)
		print "genome %s: Added !" % genome
	else:
		print "no Genome found ! O.o Weirrrrddddd"	

# raw_input("ready ?")

n = 0
running = {}
process = {}
name = {}
path = {}
UID = 0
log = open('gene_discovery.log','w')
while cmdToRun:
	try:
		while n < nb_CPU:
			UID += 1
			newcmd = cmdToRun.pop()
			print " ".join(newcmd)
			if "augustus" in newcmd:
				path[UID] = newcmd[3].split('=')[1]
			elif "snap" in newcmd:
				path[UID] = newcmd[4]
			P = Popen(newcmd, stdout=log, stderr=log)
			print "Launched %s/%s"%(UID, len(cmdToRun))
			running[UID] = 1
			process[UID] = P
			name[UID] = " ".join(newcmd)
			n += 1
		for UID in running.keys():
			if process[UID].poll() == 0:
				del running[UID]
				print "FINISHED" + name[UID]
				n -= 1
		time.sleep(1)
	except KeyboardInterrupt:
		for k in running.keys():
			os.remove(path[UID])
		print "cleaned running file"
		print "exiting"

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
if not os.path.isdir('SLURM'):
	os.mkdir('SLURM')

jobID = 0 
for cmd in cmdToRun:
	f = open('./SLURM/CDS_predict_%s.sh'%jobID,'w')
	f.write('#!/bin/sh\n')
	f.write(' '.join(cmd))
	f.close()
	jobID += 1

f = open('SLURM_CDS_prediction.sh','w')
f.write('''#!/bin/bash 
# 
#SBATCH -J tophat # A single job name for the array 
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # All cores on one machine 
#SBATCH -p serial_requeue # Partition 
#SBATCH --mem 4000 # Memory request 
#SBATCH -t 0-24:00 # 2 hours (D-HH:MM) 
#SBATCH -o log_CDS_discovery_%A_%a.out # Standard output 
#SBATCH -e log_CDS_discovery_%A_%a.err # Standard error

module load centos6/augustus-3.0
module load centos6/snap-2013-11-29''')

f.write('\n\nsh %s/SLURM/CDS_predict_"${SLURM_ARRAY_TASK_ID}".sh' % os.path.abspath('.'))
f.close()
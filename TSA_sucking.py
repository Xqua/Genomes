#!/usr/bin/env python 

import pandas as pd
import os, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the genomes ID in the collumn named 'Accession number'")
parser.add_option("-t", "--file_type", dest="filetype", default="fasta",
                  help="Type of file to be downloaded, fasta of genbank (valid option: 'fasta', 'gb' )")
(options, args) = parser.parse_args()

pathcsv = options.pathcsv
filetype = options.filetype

if pathcsv == "None":
	print "List of genome to dowload must be provided.\n -h for more information"
	sys.exit(1)

baseURL = "ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/"
datatofetch = pd.read_csv(pathcsv)
tmp = datatofetch['Accession number'].values
gIDtofetch = []
for i in tmp:
	letters = i[0:4]
	nb = i[-1]
	if filetype == 'fasta':
		TSAID = "tsa.%s.%s.fsa_nt.gz" % (letters, nb)
	elif filetype == 'gb':
		TSAID = "tsa.%s.%s.gbff.gz" % (letters, nb)
	else:
		print "Filetype not recognize, exiting"
		sys.exit(1)
	gIDtofetch.append(TSAID)

for ID in gIDtofetch:
    os.system('wget -r --no-remove-listing %s%s' % (baseURL, ID))


#!/usr/bin/env python 

import pandas as pd
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the genomes ID in the collumn named 'Accession number'")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv

if pathcsv == "None":
	print "List of genome to dowload must be provided.\n -h for more information"
	sys.exit(1)

baseURL = "ftp://ftp.ncbi.nih.gov/genomes/all/"
datatofetch = pd.read_csv(pathcsv)
gIDtofetch = datatofetch['Accession number'].values


f = open('AllGenomes')
lines = f.readlines()
genome_list = {}
for l in lines:
    if l:
        genome_list[l.strip().lower()] = l.strip()

tofetch = []
for i in gIDtofetch:
    for j in genome_list.keys():
        if i.lower() in j:
            tofetch.append(genome_list[j])

for ID in tofetch:
    os.system('wget -r --no-remove-listing %s%s' % (baseURL, ID))


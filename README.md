# SCRIPTS

## genome_sucking.py
Spider to crawl the NCBI genome FTP database.
Needs to be provided with a CSV file containing a collumn named:
"Accession number"

This collumn will be used to match the Genome ID with the NCBI database

This script will populate the folder ftp.ncbi.nih.gov/genomes/all with all the genomes.
```
Usage: genome_sucking.py [options]

Options:
  -h, --help            show this help message and exit
  -p PATHCSV, --csv_file=PATHCSV
                        [Required] Location of the csv file containing all the
                        genomes ID in the collumn named 'Accession number'
```

## Gene_Prediction_Genome_DB.py

Script that will perform CDS discovery to the genomes using different search algorithms.
Current version uses SNAP and Augustus (should be extended to using GeneID and GeneMark 
for complete discovery)

This script takes as input the location of the downloaded genomes (ex: ftp.ncbi.nih.gov/genomes/all )
It will then crawl through them and launch augustus and snap on them, and gunzip them if needed. 
You can set the number of CPU to be used for the analysis. 

Current version uses the Fly ORF/CDS HMM model.

more info on SNAP: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC421630/

more info on Augustus: http://bioinf.uni-greifswald.de/augustus/
```
Usage: Gene_Prediction_Genome_DB.py [options]

Options:
  -h, --help            show this help message and exit
  -g GENOMEFOLDER, --genome_folder=GENOMEFOLDER
                        [Required] Location of the folder containing all the
                        genomes
  -c CPU, --cpu=CPU     Number of CPU to use for the analysis. Each core is
                        used to launch one instance of SNAP or Augustus.
```

## Build_nr_DB.py

Script that crawl the genomes (--genome) and looks for the output of SNAP and Augustus.

If both files are present, it extracts the protein sequences, merge them into one fasta file.
This fasta file is then sorted by sequence lenght, and uClust is launched on it to 
create a non redundant database (nrdb).

When uClust created a nrdb for all available files, it then annotates them with the genome name.
Those files are then merged into one big nrdb file (--database) containing all the protein sequences. 

You can set the overlap threshold at which uClust will cluster two sequences together
using the parameter (--precision).

This script takes as input the location of the downloaded genomes (ex: ftp.ncbi.nih.gov/genomes/all )

more info on uClust: http://drive5.com/usearch/manual/uclust_algo.html
```
Usage: Build_nr_DB.py [options]

Options:
  -h, --help            show this help message and exit
  -g GENOMEFOLDER, --genome_folder=GENOMEFOLDER
                        [Required] Location of the folder containing all the
                        genomes
  -d DB_PATH, --database_path=DB_PATH
                        Location of the output file that will contain the non
                        redundant database
  -a, --append          Append to the database (This will not erase the DB)
  -p UCLUST_PRECISION, --precision=UCLUST_PRECISION
                        Overlap Threshold for uClust to call two sequence
                        identical (uClust will keep the longer sequence of the
                        two)
```

# REQUIREMENT

SNAP software 
included in the repository (hard to find and not maintained anymore)

Augustus Software

http://bioinf.uni-greifswald.de/augustus/downloads/

http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz

uClust

http://www.drive5.com/uclust/downloads1_2_22q.html

Pandas

http://pandas.pydata.org/

can be installed on a linux machine with:
`sudo pip install pandas`



# MISC

for more information on gene discovery 
https://www.broadinstitute.org/annotation/genome/neurospora/GeneFinding.html

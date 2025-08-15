Volvox_scan Software: Release Version 1.0

Program Description:
-------------------
The volvox_scan DNA mutation clustering software is a standalone program with code written entirely using Python language (version 3). The program predicts likely false-negative DNA mutations and likely related mutants from a sequenced population of independently mutagenized individuals. The program utilizes variant call format (VCF) files containing unfiltered, initial variant-calls from single-individual genotyping calls. The program does not currently support VCF files for joint-genotyping calls. By default the program requires no input parameters and searches for any collection or mitxure of commpressed (zipped) and/or uncompressed VCF files in the current folder (or directory) where the program is executed or invoked. On successful program analysis completion the results will be located in a "Results" folder. The program is platform-independent and has been tested successfully on Windows, UNIX/Linux and Mac OS.  

Example Runs:
------------
python3 volvox_scan.py --snps -d ../data/both/dupl
python3 volvox_scan.py --all -q 10 -d ../data/both/dupl 
python3 volvox_scan.py --all -q 12 -d ../data/v2.1_255/ -n 5 -m 250 
python3 volvox_scan.py --all -q 20 --uniq -d ../data/v2.1_255/ -n 5 -m 250 

usage: volvox_scan.py [-h] [-d DIR] [-p PREFIX] [-q QUAL] [-m MIN]
                              [-n MAXPOP] [--snps] [--all] [--uniq]

Detection of Clusters of Likely Related Mutants in a Sequenced Population

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --dir DIR     Path For Directory Containing VCF Files
  -p PREFIX, --prefix PREFIX
                        Prefix For Output File Names
  -q QUAL, --qual QUAL  Minimum Variant Quality for Mutation Spectrum, Heterozygosity and Common Variant Filtering Analyses
  -m MIN, --min MIN     Minimum No. of Shared Variants to Report a Cluster
  -n MAXPOP, --maxpop MAXPOP
                        Maximum No. of Individuals in a Cluster
  --snps                Use SNPs Only in Cluster Analysis
  --all                 Use All Variants in Cluster Analysis
  --uniq                Perform Common Variant Subtraction

Note:
====
The parameter -q specifies the minimum variant quality-score for mutations used in the  mutation spectrum and percent heterozygosity estimations for the shared DNA variants only. It is also used in filtering unique variants from the common variant subtraction option (--uniq).

#Example Runs:
-------------

RUN 1:

	Command: 
		python3 volvox_scan.py

	Explanation:
		-Search current folder for VCF files
		-Perform clustering using default parameters 
		-i.e. Default q=10, n=5, m=250, --snps
		-Retain clusters containing at most five mutant individuals sharing at least 250 mutations 
		-Use Q10 mutations or higher for the mutation spectrum and heterozygosity calculations.

RUN 2:

	Command: 
		python3 volvox_scan.py --snps -d ../data/both/dupl

	Explanation:
		-Search subfolder ../data/both/dupl in the current directory for VCF files
		-Cluster using only SNPs (i.e. exclude INDELs in the analysis)
		-Default q=10, n=5, m=250
		
RUN 3:

	Command: 
		python3 volvox_scan.py --all -q 12 -d data/both/dupl >& xn.txt

	Explanation:
		-Search subfolder data/both/dupl in the current directory for VCF files
		-Cluster using all mutations (i.e both SNPs and INDELs) in the analysis
		-Set variant quality-score threshold (-q) to 12
		-Default n=5, m=250
		-Redirect console screen output to a text file called xn.txt 

RUN 4:

	Command: 
		python3 volvox_scan.py --all -q 20 -d ../data/v2.1_255/ -n 4 -m 300 

	Explanation:
		-Search subfolder ../data/v2.1_255  in the current directory for VCF files
		-Cluster using all mutations (i.e both SNPs and INDELs) in the analysis
		-Set variant quality-score threshold (-q) to 20
		-Use Q20 mutations or higher for the mutation spectrum and heterozygosity calculations.
		-Retain clusters of at most four mutant individuals sharing at least 300 mutations 

RUN 5:

	Command: 
		python3 volvox_scan.py --all -q 15 --uniq -d .  -n 4 -m 500 

	Explanation:
		-Search current directory for VCF files
		-Cluster using all mutations (i.e both SNPs and INDELs) in the analysis
		-Perform common variant subtraction
		-Set variant quality-score threshold (-q) to 15
		-Use Q15 mutations or higher for the mutation spectrum and heterozygosity calculations.
		-Retain clusters of at most four mutant individuals sharing at least 500 mutations 
		-Retain unique mutations with Q15 or higher quality after common variants subtraction.



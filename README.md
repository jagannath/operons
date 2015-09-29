# operons
Project understanding order of genes in bacterial operons

Genomes downloaded from : ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/
File downloaded: all.ptt.tar.gz
Date downloaded: 25th September 2015


Information about the genome
Total number of bacterial species = 2,774
Total number of genomes (i.e chromosomal + plasmid) = 5,220



1. parsePTT.py: This python file parses the .ptt files (containing information of the proteins in the genome). The ptt files have COG mappings as well. 
Example of the format - 
	Acaryochloris marina MBIC11017 chromosome, complete genome - 1..6503724
	6254 proteins
	Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
	1627..2319	-	230	158333234	-	AM1_0003	-	COG1051F	NUDIX hydrolase

The Synonym is the locus_tag; The PID is the protein ID; This can be mapped again to the COG database if needed.
The file outputs (1) a dictionary of COG: [locus_tag1, locus_tag2, ...] and (2) a dictionary of (COG1,COG2): Count

2. graphCOG.py: This python file parses the COGpair dictionary and calculates the FGOC and dirFGOC scores. It makes appropriately formatted tab delimited file for use in cytoscape

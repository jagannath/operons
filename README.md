# operons
Project understanding order of genes in bacterial operons

Genomes downloaded from : ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/
**
File downloaded: all.ptt.tar.gz
Date downloaded: 25th September 2015
**
**
File downloaded: all.fna.tar.gz
Date downloaded: 8th November 2015
**


Information about the genome
Total number of bacterial species = 2,774
Total number of genomes (i.e chromosomal + plasmid) = 5,220



1. parsePTT.py: This python file parses the .ptt files (containing information of the proteins in the genome). The ptt files have COG mappings as well. Relevant dictionaries with different mappings are created such as COG:[locusTags],organism:[COG pairs] etc. 

The Synonym is the locus_tag; The PID is the protein ID; This can be mapped again to the COG database if needed.
The file outputs (1) a dictionary of COG: [locus_tag1, locus_tag2, ...] and (2) a dictionary of (COG1,COG2): Count

2. graphCOG.py: This python file parses the COGpair dictionary and calculates the FGOC and dirFGOC scores. It makes appropriately formatted tab delimited file for use in cytoscape. The notebook file also contains computation to better understand the graph. (a) Histogram of number of edges targeting and emanating from each node and also the FGOC scores. Some information of the most connected COG groups in the graph (depending on the cutoffs) is also ascertained. 

3. operonStats: This notebook contains the script for understanding the membership pattern of gene pairs within and between operons. Taking E.coli and B.subtilis (strains) as model organisms with known and annotated operon membership, plot of gene pair membership within or between operon as a function of the FGOC score is calculated. I also try to understand the correlation between intergenetic distances and the fgoc scores. 
